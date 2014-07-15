from KLog import *
lgx =logging.getLogger("KalmanFilter")
lgx.setLevel(logging.WARN)
lgx.addHandler(ch)
debug = Debug.info.value