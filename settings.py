# Avi Persin, Revision 2016-01-06

import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

BASE_PATH = BASE_DIR + '/gistemp4.0'

TMP_DIR = BASE_PATH + '/tmp/'

PROGRESS_DIR = TMP_DIR + 'progress/'

SOURCES_DIR = BASE_PATH + '/config/'

INPUT_DIR = TMP_DIR + 'input/'

LOG_DIR = TMP_DIR + 'log/'

RESULT_DIR = TMP_DIR + 'result/'

WORK_DIR = TMP_DIR + 'work/'
