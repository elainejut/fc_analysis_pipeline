# FC Analysis Pipeline 

This repository contains the code for the analysis pipeline that compares pre- and post-psychedelic resting-state EEG FC.

## Usage

1) Open "main.mlx" --> define all general variables at the top of the file
2) Main parameters to think about:

	a) FC_METHOD --> options: ['icoh', 'amplcorr', 'mi', 'dtf', 'pdc']
	
	b) RUN_FROM_BEGINNING --> if FC analysis has already been run on some data, this can be set to false to continue running from the saved data rather than restarting
	
	c) REREF_MODE --> re-referencing mode: "car"=common average referencing; "none"=no re-referencing is applied
	
	d) DOWNSAMPLE_BOOL --> "true": downsample the data to 250 Hz; "false": keep the data at 500 Hz
3) Define "SAVE_DIR" in the beginning, and all data and figures will be saved into this directory
4) MODEL_ORDERS_MAT --> right now the model orders are set to the predefined data - DD1BL vs DD2150; if switching to different dataset, either set MVAR model order as a different value or find the individual model orders per data sample
5) After defining the parameters at the top of the file, run "main.mlx" 

*For any questions, feel free to contact ej24@ic.ac.uk
	