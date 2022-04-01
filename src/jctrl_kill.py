import os
for id in range(136008, 136082):
		os.popen('jctrl kill '+str(id))
		os.popen('sleep 1')
