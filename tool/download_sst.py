# download_sst.py
# make sure to update paths in toolconfig.py first

import datetime
import requests
import subprocess
import os
import time
import json
import sys

from bs4 import BeautifulSoup
import toolconfig

CURL_PATH = str(toolconfig.CURL_PATH)
SST_NCEI_PATH = toolconfig.SST_NCEI_PATH
BASE_URL = "https://www.ncei.noaa.gov/data/sea-surface-temperature-extended-reconstructed/v6/access/"
MANIFEST_FILE = str(SST_NCEI_PATH / "download_manifest.json")
TIME_WINDOW = 48 * 3600  # 48 hours in seconds


# ---------------------------------------------------------

def load_manifest():
    try:
        with open(MANIFEST_FILE, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        return {}

def save_manifest(manifest):
    with open(MANIFEST_FILE, 'w') as f:
        json.dump(manifest, f, indent=2)


# ---------------------- HTML PARSER ----------------------

def get_server_listing():
    """Pull files from NOAA html directory listing"""
    print("Fetching server listing...")
    r = requests.get(BASE_URL, timeout=120)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, 'html.parser')
    rows = soup.find_all("tr")

    files = {}

    for row in rows:
        cols = row.find_all("td")
        if len(cols) < 3:
            continue

        link = cols[0].find("a")
        if not link:
            continue

        filename = link.text.strip()
        if not filename.endswith(".nc"):
            continue

        timestr = cols[1].text.strip()
        sizestr = cols[2].text.strip()

        try:
            size = int(sizestr)
            server_time = datetime.datetime.strptime(timestr, "%Y-%m-%d %H:%M").replace(tzinfo=datetime.timezone.utc)
            files[filename] = {
                "url": BASE_URL + filename,
                "size": size,
                "mtime": int(server_time.timestamp())
            }
        except:
            continue

    return files


# ---------------------- CHECKS ----------------------

def is_time_match(local_mtime, server_mtime):
    return abs(local_mtime - server_mtime) <= TIME_WINDOW


def is_valid_file(filename, server_info, manifest):
    if not os.path.exists(filename):
        return False

    local_size = os.path.getsize(filename)
    local_time = int(os.path.getmtime(filename))

    if local_size != server_info["size"]:
        return False

    if not is_time_match(local_time, server_info["mtime"]):
        return False

    if filename in manifest:
        if manifest[filename]["size"] != local_size:
            return False
        if not is_time_match(manifest[filename]["mtime"], server_info["mtime"]):
            return False

    return True


# ---------------------- DOWNLOAD ----------------------

def download_file(url, filename, expected_size, expected_mtime):
    """
    Download with curl. After verifying size, set file mtime to expected_mtime (server mtime).
    Returns True on success, False otherwise.
    """
    for _ in range(3):
        try:
            subprocess.run(
                [CURL_PATH, "--max-time", "120", "-fLo", filename, url],
                check=True
            )

            if os.path.exists(filename) and os.path.getsize(filename) == expected_size:
                # set atime and mtime to server mtime (seconds since epoch)
                try:
                    os.utime(filename, (expected_mtime, expected_mtime))
                except Exception as e:
                    print(f"Warning: failed to set mtime for {filename}: {e}")
                return True

            # size mismatched -> remove and retry
            if os.path.exists(filename):
                os.remove(filename)

        except Exception:
            if os.path.exists(filename):
                os.remove(filename)
            time.sleep(15)

    return False


# ---------------------- MAIN ----------------------

def main():
    server_files = get_server_listing()
    if not server_files:
        print("No files discovered from server.")
        sys.exit(1)

    manifest = load_manifest()

    if not manifest:
        print("Creating new manifest from existing files (if present)...")

        for remotefname, info in server_files.items():
            fname = str(SST_NCEI_PATH / remotefname)
            if os.path.exists(fname) and os.path.getsize(fname) == info['size']:
                # if existing file matches size but mtime may differ -> keep local mtime only if within TIME_WINDOW
                local_mtime = int(os.path.getmtime(fname))
                if is_time_match(local_mtime, info['mtime']):
                    manifest[fname] = {
                        "path": os.path.abspath(fname),
                        "size": os.path.getsize(fname),
                        "mtime": int(local_mtime)
                    }
                else:
                    # prefer server mtime only after a fresh download; treat as mismatch for now
                    continue

        save_manifest(manifest)

    failed_queue = []

    for remotefname, info in sorted(server_files.items()):
        fname = str(SST_NCEI_PATH / remotefname)
        url = info['url']

        if is_valid_file(fname, info, manifest):
            print(f"✔ Skipping valid: {fname}")
            continue

        if os.path.exists(fname):
            print(f"✖ Deleting mismatched: {fname}")
            os.remove(fname)

        print(f"⬇ Downloading {fname}")
        success = download_file(url, fname, info['size'], info['mtime'])

        if success:
            # store server mtime in manifest (not ephemeral local mtime)
            manifest[fname] = {
                "path": os.path.abspath(fname),
                "size": info["size"],
                "mtime": info["mtime"]
            }
            save_manifest(manifest)
        else:
            print(f"⚠ Queuing for retry: {fname}")
            failed_queue.append((url, fname, info["size"], info["mtime"]))

    # ---------------------- RETRY QUEUE ----------------------

    if failed_queue:
        print(f"\n🔁 Retrying {len(failed_queue)} failures up to 5 rounds\n")

        for attempt in range(5):
            print(f"Retry round {attempt+1}")

            still_failed = []

            for url, fname, size, mtime in failed_queue:
                if os.path.exists(fname):
                    # if somehow present, check size and set mtime if needed
                    if os.path.getsize(fname) == size:
                        try:
                            os.utime(fname, (mtime, mtime))
                        except Exception as e:
                            print(f"Warning: failed to set mtime for {fname} during retry: {e}")
                        manifest[fname] = {
                            "path": os.path.abspath(fname),
                            "size": size,
                            "mtime": mtime
                        }
                        save_manifest(manifest)
                        continue
                    else:
                        os.remove(fname)

                print(f"Retrying: {fname}")
                if download_file(url, fname, size, mtime):
                    manifest[fname] = {
                        "path": os.path.abspath(fname),
                        "size": size,
                        "mtime": mtime
                    }
                    save_manifest(manifest)
                else:
                    still_failed.append((url, fname, size, mtime))

            if not still_failed:
                break

            failed_queue = still_failed
            time.sleep(20)

    if failed_queue:
        print("\n❌ Some files failed permanently:")
        for _, fname, _, _ in failed_queue:
            print("   -", fname)
        sys.exit(1)

    print("\n✅ All files verified and synced")
    sys.exit(0)


if __name__ == "__main__":
    main()
