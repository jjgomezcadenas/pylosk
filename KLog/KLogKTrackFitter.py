from KLog import *

lgx =logging.getLogger("KTrackFitter")
lgx.setLevel(logging.DEBUG)
lgx.addHandler(ch)
debug = Debug.info.value
