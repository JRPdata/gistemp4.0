"""
Equal-area grid geometry, translated as literally as possible from:
  - NTRP2EA0 / NTRP2EA subroutines in MaskRegrid.f
  - def_ea_grid subroutine in rearrange_ERSST.f

All arrays below use 1-based indexing (index 0 is left unused / zero) so that
the loop bodies mirror the Fortran source line-for-line and are easy to
audit against the .f files.
"""
import math
import numpy as np

TWOPI = 6.283185307179586477

IMA = 180   # input (ERSST) grid: longitude boxes
JMA = 89    # input (ERSST) grid: latitude boxes
OFFIA = 89.5
DIVJA = 90.0
SKIB = 9999.0
JMB = 80    # output (equal-area) grid: number of latitude bands


def nint(x):
    """Fortran NINT: round to nearest integer, ties away from zero."""
    if x >= 0:
        return int(math.floor(x + 0.5))
    else:
        return int(math.ceil(x - 0.5))


class EAGrid:
    """
    Holds everything needed to regrid a (89,180) [lat,lon] ERSST array onto
    the 8000-box Sergei equal-area grid, plus the box ordering/boundaries
    used by the final SBBX file.
    """

    def __init__(self):
        self._build_ntrp2ea0()
        self._build_edges()
        self._build_def_ea_grid()

    # ------------------------------------------------------------------
    # NTRP2EA0 (MaskRegrid.f) -- one-time grid partition setup
    # ------------------------------------------------------------------
    def _build_ntrp2ea0(self):
        IMA, JMA = self.IMA, self.JMA = 180, 89

        if not (1 <= IMA <= 720) or not (1 <= JMA <= 361):
            raise ValueError("invalid grid dimensions")

        # ---- I direction partition (longitude), per izone = 1..4 -------
        # IMIN/IMAX/FMIN/FMAX indexed [IB(1..160), izone(1..4)]
        IMIN = np.zeros((161, 5), dtype=np.int64)
        IMAX = np.zeros((161, 5), dtype=np.int64)
        FMIN = np.zeros((161, 5), dtype=np.float64)
        FMAX = np.zeros((161, 5), dtype=np.float64)

        DIA = 360.0 / IMA
        for izone in range(1, 5):
            IMB = 40 * izone
            DIB = 360.0 / IMB
            IA = 1
            RIA = (IA + OFFIA) * DIA - 360
            IB = IMB
            for IBP1 in range(1, IMB + 1):
                RIB = (IBP1 - 1) * DIB
                while True:
                    diff = RIA - RIB
                    if diff < 0:
                        IA += 1
                        RIA += DIA
                        continue
                    elif diff == 0:
                        IMAX[IB, izone] = IA
                        FMAX[IB, izone] = 0.0
                        IA += 1
                        RIA += DIA
                        IMIN[IBP1, izone] = IA
                        FMIN[IBP1, izone] = 0.0
                        break
                    else:
                        IMAX[IB, izone] = IA
                        FMAX[IB, izone] = (RIA - RIB) / DIA
                        IMIN[IBP1, izone] = IA
                        FMIN[IBP1, izone] = 1.0 - FMAX[IB, izone]
                        break
                IB = IBP1
            IMAX[IMB, izone] += IMA

        self.IMIN, self.IMAX, self.FMIN, self.FMAX = IMIN, IMAX, FMIN, FMAX

        # ---- J direction partition (latitude), input grid A -----------
        SINA = np.zeros(JMA + 1, dtype=np.float64)
        OFFJA = (DIVJA - JMA) / 2.0
        DJA = 0.5 * TWOPI / DIVJA
        for JA in range(1, JMA):  # 1 .. JMA-1
            RJA = (JA + OFFJA) * DJA - 0.25 * TWOPI
            SINA[JA] = math.sin(RJA)
        SINA[0] = -1.0
        SINA[JMA] = 1.0
        self.SINA = SINA

        # ---- J direction, output equal-area grid B ---------------------
        SINB = np.zeros(81, dtype=np.float64)
        SINB[0] = -1.0
        sband = SINB[0]
        for izone in range(1, 5):
            dsband = 0.1 * izone
            for jzs in range(1, 11):
                jb = jzs + (izone - 1) * 10
                SINB[jb] = sband + 0.1 * jzs * dsband
            sband += dsband
        SINB[40] = 0.0
        for jb in range(41, 81):
            SINB[jb] = -SINB[80 - jb]
        self.SINB = SINB

        # ---- nbefor: cumulative box count south of each band -----------
        nbefor = np.zeros(82, dtype=np.int64)
        jb = 1
        nbefor[1] = 0
        for izone in range(1, 9):
            iz = izone if izone <= 4 else 9 - izone
            for jzs in range(1, 11):
                nbefor[jb + 1] = nbefor[jb] + 40 * iz
                jb += 1
        self.nbefor = nbefor

        # ---- JMIN/JMAX/GMIN/GMAX ----------------------------------------
        JMIN = np.zeros(82, dtype=np.int64)
        JMAX = np.zeros(81, dtype=np.int64)
        GMIN = np.zeros(82, dtype=np.float64)
        GMAX = np.zeros(81, dtype=np.float64)
        JMIN[1] = 1
        GMIN[1] = 0.0
        JA = 1
        for JB in range(1, JMB):  # 1 .. JMB-1
            while True:
                diff = SINA[JA] - SINB[JB]
                if diff < 0:
                    JA += 1
                    continue
                elif diff == 0:
                    JMAX[JB] = JA
                    GMAX[JB] = 0.0
                    JA += 1
                    JMIN[JB + 1] = JA
                    GMIN[JB + 1] = 0.0
                    break
                else:
                    JMAX[JB] = JA
                    GMAX[JB] = SINA[JA] - SINB[JB]
                    JMIN[JB + 1] = JA
                    GMIN[JB + 1] = SINB[JB] - SINA[JA - 1]
                    break
        JMAX[JMB] = JMA
        GMAX[JMB] = 0.0

        self.JMIN, self.JMAX, self.GMIN, self.GMAX = JMIN, JMAX, GMIN, GMAX

    # ------------------------------------------------------------------
    # Precompute (box, ia, ja, coef) contribution edges, and per-box
    # total_wt (== sum of F*G over the box's edges).
    #
    # total_wt is accumulated in float32, in the exact same order as the
    # Fortran DO 510 loops, because `total_wt` is declared REAL*4 in
    # MaskRegrid.f (rounded after every addition) while WEIGHT/VALUE are
    # REAL*8. This only affects the WEIGHT > total_wt*frac gate, never the
    # output SST value itself (which comes from the REAL*8 VALUE/WEIGHT).
    # ------------------------------------------------------------------
    def _build_edges(self):
        IMA = self.IMA
        IMIN, IMAX, FMIN, FMAX = self.IMIN, self.IMAX, self.FMIN, self.FMAX
        JMIN, JMAX, GMIN, GMAX = self.JMIN, self.JMAX, self.GMIN, self.GMAX
        SINA, nbefor = self.SINA, self.nbefor

        edge_box = []
        edge_ia = []
        edge_ja = []
        edge_coef = []
        total_wt = np.zeros(8000, dtype=np.float32)

        for JB in range(1, JMB + 1):
            izone = (JB + 9) // 10
            if izone > 4:
                izone = 9 - izone
            imb = izone * 40
            JAMIN = JMIN[JB]
            JAMAX = JMAX[JB]
            for IB in range(1, imb + 1):
                IJB = IB + nbefor[JB]
                box_idx = IJB - 1
                IAMIN = IMIN[IB, izone]
                IAMAX = IMAX[IB, izone]
                box_total_wt = np.float32(0.0)
                for JA in range(JAMIN, JAMAX + 1):
                    G = SINA[JA] - SINA[JA - 1]
                    if JA == JAMIN:
                        G -= GMIN[JB]
                    if JA == JAMAX:
                        G -= GMAX[JB]
                    for IAREV in range(IAMIN, IAMAX + 1):
                        IA = 1 + (IAREV - 1) % IMA
                        F = 1.0
                        if IAREV == IAMIN:
                            F -= FMIN[IB, izone]
                        if IAREV == IAMAX:
                            F -= FMAX[IB, izone]
                        coef = F * G
                        edge_box.append(box_idx)
                        edge_ia.append(IA - 1)
                        edge_ja.append(JA - 1)
                        edge_coef.append(coef)
                        box_total_wt = np.float32(box_total_wt + np.float32(coef))
                total_wt[box_idx] = box_total_wt

        self.edge_box = np.asarray(edge_box, dtype=np.int64)
        self.edge_ia = np.asarray(edge_ia, dtype=np.int64)
        self.edge_ja = np.asarray(edge_ja, dtype=np.int64)
        self.edge_coef = np.asarray(edge_coef, dtype=np.float64)
        self.total_wt = total_wt

    # ------------------------------------------------------------------
    # NTRP2EA entry point: interpolate one (89,180) [lat,lon] map onto the
    # 8000-box equal-area grid ("natural" south->north, west->east order,
    # 0-based here; this is the order MaskRegrid.f writes to ERdSST_monthly).
    # ------------------------------------------------------------------
    def regrid(self, data, wta, frac=0.50, skip=9999.0):
        """
        data, wta: numpy arrays shaped (89,180) == [lat index 0..88, lon index 0..179]
        Returns: float32 array of length 8000 (natural south->north/west->east order).
        """
        if frac > 0.99:
            frac = 0.99
        if frac < 0.0:
            frac = 0.0

        data_edge = data[self.edge_ja, self.edge_ia]
        wta_edge = wta[self.edge_ja, self.edge_ia]
        valid = data_edge != np.float32(skip)

        coef = self.edge_coef
        w_c = np.where(valid, coef * wta_edge.astype(np.float64), 0.0)
        v_c = np.where(valid, w_c * data_edge.astype(np.float64), 0.0)

        WEIGHT = np.bincount(self.edge_box, weights=w_c, minlength=8000)
        VALUE = np.bincount(self.edge_box, weights=v_c, minlength=8000)

        thresh = self.total_wt.astype(np.float64) * frac
        good = WEIGHT > thresh
        safe_weight = np.where(WEIGHT == 0, 1.0, WEIGHT)
        ratio = VALUE / safe_weight

        B = np.full(8000, skip, dtype=np.float32)
        B[good] = ratio[good].astype(np.float32)
        return B

    # ------------------------------------------------------------------
    # def_ea_grid (rearrange_ERSST.f): box boundaries (in "natural" order)
    # and the natural-order -> Sergei-canonical-order permutation.
    # ------------------------------------------------------------------
    def _build_def_ea_grid(self):
        xbypi = 9000.0 / math.asin(1.0)  # 100 * 180/pi

        nbefor = np.zeros(82, dtype=np.int64)
        LatSb = np.zeros(82, dtype=np.int64)
        j = 1
        nbefor[1] = 0
        sband = -1.0
        for iband in range(1, 9):
            iz = iband if iband <= 4 else 9 - iband
            dsband = 0.1 * iz
            for jzs in range(1, 11):
                nbefor[j + 1] = nbefor[j] + 40 * iz
                LatSb[j] = nint(xbypi * math.asin(sband + 0.1 * (jzs - 1) * dsband))
                j += 1
            sband += dsband
        LatSb[81] = -LatSb[1]

        # sanity: this nbefor must equal the one built in _build_ntrp2ea0
        assert np.array_equal(nbefor, self.nbefor), "nbefor mismatch between NTRP2EA0 and def_ea_grid"

        # latlon[k] = [lts, ltn, lnw, lne], k = 1..8000 (natural order)
        latlon = np.zeros((8001, 4), dtype=np.int64)
        for jband in range(1, 81):
            latsouth = LatSb[jband]
            latnorth = LatSb[jband + 1]
            im = nbefor[jband + 1] - nbefor[jband]
            for i in range(1, im + 1):
                lonwest = -18000 + nint((i - 1) * 36000.0 / im)
                loneast = -18000 + nint(i * 36000.0 / im)
                idx = nbefor[jband] + i
                latlon[idx] = [latsouth, latnorth, lonwest, loneast]

        # ij_SN_WE(ij_sb) = SN_WE (natural) index for canonical position ij_sb
        ij_sn_we = np.zeros(8001, dtype=np.int64)
        ij = 0
        for jband in range(1, 9):
            jup = jband * 10 + 1
            iz = jband if jband <= 4 else 9 - jband
            nbefore = 8000 - nbefor[jup]
            for j in range(1, 11):
                for nb in range(1, 4 * iz + 1):
                    nbeforj = nbefore + (nb - 1) * 100
                    for i in range(1, 11):
                        ij += 1
                        idx = nbeforj + (j - 1) * 10 + i
                        ij_sn_we[idx] = ij

        self.latlon = latlon         # 1-based, index 1..8000 -> [lts,ltn,lnw,lne]
        self.ij_sn_we = ij_sn_we     # 1-based, index 1..8000 (canonical n) -> natural index
