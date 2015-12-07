
import sys 

''' Number of lines after which to update progress counter '''
PROGRESS_DELAY = 1000



############################################################
# Progress counter
############################################################
def update_progress(progress, operation='Progress',  max_iter=None, skip_lines=PROGRESS_DELAY):
    '''
	update_progress(): Displays or updates a console progress bar
	
    USAGE: update_progress(progress, operation='Progress',  max_iters=None, progress_delay=PROGRESS_DELAY)

    @arg progress: Accepts a float between 0 and 1. Any int will be converted to a float.
    A value under 0 represents a 'halt'.
    A value at 1 or bigger represents 100%

    @arg operation: What to show before progress bar. 
	
    @arg max_iter: Alternatively if max_iter is specified - then progress is interpreted as current_iteration. 

    @arg skip_lines: Number of lines to skip before updating. 
    '''
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    update_progress_flag = False 
    if (max_iter is None):
        update_progress_flag = True
        if isinstance(progress, int):
			progress = float(progress)
        if not isinstance(progress, float):
			progress = 0
			status = "error: progress var must be float\r\n"
        if progress < 0:
			progress = 0
			status = "Halt...\r\n"
        if progress >= 1:
			progress = 1
			status = "Done...\r\n"
    elif (max_iter > skip_lines) and ( progress % skip_lines is 0): 
        update_progress_flag = True
        progress = float(progress) / max_iter
	## do not update each iteration 
    if update_progress_flag is True:
		block = int(round(barLength*progress))
		text = "\r{0: <40}: [{1}] {2}% {3}".format(operation, "#"*block + "-"*(barLength-block), progress*100, status)
		sys.stderr.write(text)
		sys.stderr.flush()
