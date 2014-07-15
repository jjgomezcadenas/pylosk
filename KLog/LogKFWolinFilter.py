from KLog import *
#create logger
lgx =logging.getLogger("KFWolinFilter")
lgx.setLevel(logging.INFO)
lgx.addHandler(ch)
debug = Debug.info.value