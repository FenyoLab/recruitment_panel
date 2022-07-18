macro "process_nd2_and_drawROIs"
{
	run("Set Measurements...", "area mean min integrated redirect=None decimal=3");
	
	dir = getDirectory("Choose the directory containing the subfolders with the nd2 files.");
	//load subfolders of the chosen directory
	subfolders=listDirs(dir);
	
	run("ROI Manager...");
	//print(subfolders);
	for (i=0; i<subfolders.length; i++)
	{// get file list from each subfolder
		
		list = getFileList(dir+'/'+subfolders[i]); 
		for (j=0; j<list.length; j++)
		{
			if (matches(list[j], ".+\.nd2"))
			{//for each nd2 file, process the file
				
				if(indexOf(list[j], 'static') == -1)
				{ 	
					
					//print(list[j]);
		
					//open nd2 file, save as tif and max z-project
					nd2_file = dir+'/'+subfolders[i]+'/'+list[j];
					tif_file=replace(nd2_file,'\.nd2$','.tif');
					lines_file=replace(nd2_file, '\.nd2$','_ROI.tif');
					roi_zip = replace(lines_file, '\.tif$','.zip');
					max_file_name=replace(list[j],'\.nd2$','.tif');
					max_file_name = "MAX_"+max_file_name;
					lines_file_name=replace(list[j],'\.nd2$','_ROI.tif');
					
					if(File.exists(roi_zip))
					{
						open(roi_zip);
						roiManager("Show All");
						open(dir+'/'+subfolders[i]+'/'+max_file_name);
					}
					else 
					{
						run("Bio-Formats Importer", "open=["+nd2_file+"] color_mode=Grayscale concatenate_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
						//waitForUser("");
						saveAs("Tiff", tif_file);
						run("Z Project...", "projection=[Max Intensity]");
						
						selectWindow(max_file_name);
						run("Measure");
						max_px=getResult("Max", 0);
						max_px=max_px+1; //should be way below 65535 or this will generate an error below...
						selectWindow("Results"); 
	         			run("Close" );
						//print(max_px);
			
						//open the ROI lines file - use BioImporter since it is 1-bit, and place selection on the max-file; then save the max-file
						run("Bio-Formats Importer", "open=["+lines_file+"] color_mode=Grayscale concatenate_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
						lines_file_name=getTitle();
						//waitForUser(lines_file_name);
						
						run("Duplicate...", "title=temp.tif");
						setThreshold(0, 0);
						run("Convert to Mask");
						run("Max...", "value=1");
						imageCalculator("Multiply",max_file_name,"temp.tif");
						
						
						selectWindow(lines_file_name);
						run("16-bit");
						run("Multiply...", "value="+max_px);
						imageCalculator("Add", max_file_name,lines_file_name);
	
						saveAs("Tiff", dir+'/'+subfolders[i]+'/'+max_file_name);
						close("\\Others");
						
					}
					selectWindow(max_file_name);
					//roiManager("Show All");
					Dialog.createNonBlocking("Select Nuclei");
					Dialog.addMessage("Draw a polygon/rectangle around each nucleus, and Add to ROI Manager (t).\nClick Okay when finished to save results. Click cancel to stop macro.");
					Dialog.show();
					//waitForUser("Select Nuclei" ,"Draw a polygon/rectangle around each nucleus, and Add to ROI Manager (t).\nOkay to continue/Cancel to stop.");
					roiManager("Show All");
					roiManager("Save", roi_zip);
					roiManager("deselect");
					roiManager("delete")
					
					close(max_file_name);
					
				}
			}
		}
	}
	selectWindow("ROI Manager");
	run("Close");
}

macro "clip_and_align_nuclei" 
{
	ce_saturate=0.1;
	filter_r=2;
	
	
	
	// goes through each tif movie and each ROI zip file
	// clips out ROI from tif movie (all time frames)
	// aligns the nuclei using plugin
	// saves result (_aligned.tif)
	
	dir = getDirectory("Choose the directory containing the subfolders with the ROI (zip) files.");
	
	run("Text Window...", "name=[Progress info] width=100 height=2");
	
	setBatchMode(true);
	
	// load subfolders of the chosen directory
	subfolders=listDirs(dir);

	for (i=0; i<subfolders.length; i++)
	{// get file list from each subfolder
		
		// create dir cropped_nuclei if does not already exist
		new_dir = dir+'/'+subfolders[i]+'/cropped_nuclei';
		if(!File.isDirectory(new_dir))
		{
			File.makeDirectory(new_dir);
		}
		
		list = getFileList(dir+'/'+subfolders[i]); 
		for (j=0; j<list.length; j++)
		{
			if (matches(list[j], ".+\_ROI.zip"))
			{// for each zip file, load tif movie also
			 // then clip out nuclei stacks, save and align

				 tif_file=replace(list[j],'\_ROI.zip$','.tif');
				 max_file='MAX_'+tif_file;
				 
				 roiManager("Open", dir+'/'+subfolders[i]+'/'+list[j]);
				 n_ROIs = roiManager("count");

				 open(dir+'/'+subfolders[i]+'/'+max_file);
				 selectImage(max_file);
				 run("Measure");
				 max_px=getResult("Max", 0);
				 max_px=max_px-1; //to set threshold
				 for(k=0; k<n_ROIs; k++)
				 {
				 	roiManager("Select", k);
				 	roi_name=Roi.getName();
				 	run("Duplicate...", " ");
				 	setThreshold(0, max_px);
				 	run("Convert to Mask");

				 	file_0=replace(list[j], '\_ROI.zip$', '-'+roi_name+'_mask.tif');
				 	saveAs("Tiff",new_dir + "/" + file_0);
				 	close();
				 	
				 	print("[Progress info]", "\\Update:Processing folder " + subfolders[i] + ", image " + list[j] + " and ROI " + k);
				 	
				 }
				 selectWindow("Results"); 
	         	 run("Close" );
				 selectImage(max_file);
				 close();
				 
				 open(dir+'/'+subfolders[i]+'/'+tif_file);
				 selectImage(tif_file);
				 for(k=0; k<n_ROIs; k++)
				 {
				 	roiManager("Select", k);
				 	roi_name=Roi.getName();
				 	type=Roi.getType();
				 	
				 	file_1=replace(list[j], '\_ROI.zip$', '-'+roi_name+'.tif');
				 	file_2=replace(list[j], '\_ROI.zip$', '-'+roi_name+'_adj.tif');
				 	file_3=replace(list[j], '\_ROI.zip$', '-'+roi_name+'_aligned.tif');

				 	run("Duplicate...", "title="+file_1+" duplicate");
				 	if(type != "rectangle")
				 	{
				 		setBackgroundColor(0, 0, 0);
				 		selectWindow(file_1);
				 		run("Clear Outside", "stack");
				 		run("Select None");
				 	}
				 	saveAs("Tiff", new_dir+"/"+file_1);
	
				 	//duplicate again and do some filter here to help alignment program
					run("Duplicate...", "title="+file_2+" duplicate");
					selectWindow(file_2);
					run("Enhance Contrast...", "saturated="+ce_saturate);
					run("Median...", "radius="+filter_r);
					saveAs("Tiff", new_dir+"/"+file_2);

					//use filtered image for alignment, but align the actual raw image
					selectWindow(file_2);
					//waitForUser("PAUSE1");
      				run("MultiStackReg", "stack_1="+file_2+" action_1=Align stack_2="+file_1+" action_2=[Align to First Stack] transformation=[Scaled Rotation]");

					selectWindow(file_1);
      				saveAs("Tiff", new_dir+"/"+file_3);
   					close();

   					selectWindow(file_2);
   					close();

				 }
				 roiManager("deselect");
				 roiManager("delete")
				 
				 selectWindow(tif_file);
				 close();
			}
		}
	}

	setBatchMode(false);
	waitForUser("Finished!");
	
}

macro "measure_stripe"
{
	PixelMax=4095; //pixel range from 0 to 4095 for 12-bit images
	ce_saturate=0.1;
	filter_r=2;
	expand_h=2;
	fps=15;
	scr_h=screenHeight;
	scr_w=screenWidth;
	animation_on=0;
	xy_scale=0.80
	
	//for each passed nuclei, it has been aligned
	//get position of initial stripe, expand into rectangle (Width=?)
	//get mean intensity inside stripe
	//detect nuclei by thresholding
	//get mean intensity in region of nuclei outside of stripe
	//save results to table...

	//load aligned nuclei image
	//load corresponding ROI file from cropped_nuclei folder

	base_dir = getDirectory("Choose the directory containing the subfolders with aligned nuclei movies.");
	subfolders=listDirs(base_dir);
	run("Clear Results");
	for (i=0; i<subfolders.length; i++)
	{// get file list from each subfolder
		
		dir = base_dir+'/'+subfolders[i]+'/cropped_nuclei';

		subfolder_name = replace(subfolders[i], "\/","");
		getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
    	str=""+year + "_"+ month+"_"+dayOfMonth+"_"+hour+"_"+minute+"_"+second;
		all_csv_name = base_dir + '/' + subfolders[i] + '/' +subfolder_name+'_'+str+'_measurements.csv';
		//listfile_name = base_dir + '/' + subfolders[i] + '/' +subfolder_name+'_'+str+'_fail_frames.txt';

		list = getFileList(dir);
		
		//start again or continue?
		ret = getBoolean(subfolder_name+": Do you want to continue from any previous work? (Answer NO if you'd like to delete results and start over.)");
		if(!ret)
		{
			for(j=0;j<list.length;j++)
			{
				if(endsWith(list[j],".roi"))
				{
					//waitForUser("Deleting "+dir+'/'+list[j]);
					File.delete(dir+'/'+list[j]);
				}
			}
		}
		
		//first_open=1;
		al_count=0;
		
		for(j=0;j<list.length;j++)
		{
			if(endsWith(list[j],"_aligned.tif"))
			{	al_count++;
				//load aligned movie and line mask file
				
				roi_selection = replace(list[j],"\_aligned.tif$","_selection.roi");
				roi_selection2 = replace(list[j],"\_aligned.tif$","_selection_nuclei.roi");
				if(File.exists(dir+'/'+roi_selection)) 
				{
					print(list[j] + " already processed, skipping...");
				}
				else 
				{
					// Load in aligned movie and mask of 'line' across nuclei
					mask_file = replace(list[j], "\_aligned.tif$","_mask.tif");
					open(dir+'/'+mask_file);
					open(dir+'/'+list[j]);
					getDimensions(width, height, channels, slices, frames);
					run("Duplicate...", "title=frame0");
					
					//display line mask as selection on tiff movie
					selectWindow(mask_file);
					run("Invert");
					run("Create Selection");
					
					
					getSelectionBounds(x,y,w,h);
					total_h=(expand_h*2)+h;
					makeRectangle(x,y-expand_h,w,total_h);
					selectWindow(list[j]);
					run("Restore Selection");
					run("In [+]");
					run("In [+]");
					run("In [+]");
					//run("In [+]");
					getLocationAndSize(x, y, width, height);
					
					//user can choose from good/bad alignment and if good alignment, set end frame (if needed) and adjust selection (if needed)
					if(animation_on)
					{
						run("Animation Options...", "speed="+fps+" first=1 last="+frames+" loop");
						doCommand("Start Animation [\\]");
					}

					Dialog.createNonBlocking("Image "+al_count);
					Dialog.setLocation(x,y+height+25);
					
					Dialog.addMessage("Uncheck if image is a bad alignment. If alignment Ok, set fail frame and adjust selection (as needed).\nPress OK when done.");
					Dialog.addCheckbox("Good alignment", 1);
					//Dialog.addNumber("FPS", fps, 0, 2, "frames per second");
					Dialog.addNumber("Expand by", expand_h, 0,2,"pixels");
					Dialog.addNumber("Fail frame", frames, 0, 3, "");
					Dialog.show();
					retVal=Dialog.getCheckbox();
					//fps=Dialog.getNumber();
					expand_h=Dialog.getNumber();
					fail_frame=Dialog.getNumber();
					
					selectWindow(list[j]);
					if(animation_on) { doCommand("Stop Animation"); }
					
					
					roiManager("reset");
					roiManager("add");
					roiManager("select",0);
					roiManager("save selected", dir+'/'+roi_selection); //do not delete, needed to know if already processed
					/////
					
					if(retVal)
					{//good alignment
	
						selectWindow(list[j]);
						getDimensions(width, height, channels, slices, frames);
						//waitForUser("PAUSE0 "+width+" "+height+" "+frames);
						
						//make mask for the selection 
						newImage("line mask", "8-bit black", width, height, 1);
						selectWindow("line mask");
						roiManager("deselect");
						roiManager("select",0);  
						Roi.setFillColor(255);
						roiManager("fill");
						//waitForUser("1");
	
	
						//make mask for the nuclei
						newImage("nuclei mask", "8-bit black", width, height, 1);
						roiManager("reset");
						
						selectWindow("frame0");
						run("Enhance Contrast...", "saturated="+ce_saturate);
						run("Median...", "radius="+filter_r);
						setAutoThreshold("Otsu dark");
						selectWindow("frame0");
						run("Convert to Mask");
						run("Fill Holes");
	
						run("Analyze Particles...", "add");
						n_ROIs=roiManager("count");
						max_area=0;
						max_area_pos=0;
						for(k=0;k<n_ROIs;k++)
						{
							roiManager("Select", k);
							Roi.getBounds(x, y, width_, height_);
							//print(x+" "+y+" "+width_+" "+height_);
							area=width_*height_;
							if(area > max_area)
							{
								max_area=area;
								max_area_pos=k;
							}
						}
						
						selectWindow("nuclei mask");
						roiManager("deselect");
						roiManager("select", max_area_pos);
						//waitForUser("2");

						//rescale the ROI down a bit in case alignment moves cell a small amount it will still be captured without blank space
						run("Scale... ", "x="+xy_scale+" y="+xy_scale+" centered");
						roiManager("Add");
						n_ROIs=roiManager("count");
						roiManager("select", n_ROIs-1); // make sure last ROI is selected - the scaled ROI
						
						roiManager("save selected", dir+'/'+roi_selection2);  
						Roi.setFillColor(255);
						roiManager("fill");
						
						roiManager("reset");
						selectWindow("ROI Manager");
						run("Close");
						//waitForUser("3");
						
						// line:
						imageCalculator("AND create", "nuclei mask","line mask");
						nuclei_line_mask_title=getTitle();
						//rename("nuclei line mask"); //02_10_21 SK
						//waitForUser(nuclei_line_mask_title);
						
						//background:
						imageCalculator("Difference", "nuclei mask",nuclei_line_mask_title); //"nuclei line mask"); //02_10_21 SK
						background_mask_title=getTitle();
						//rename("background mask"); //02_10_21 SK
						//waitForUser(background_mask_title);
						
						//background
						selectWindow(list[j]);
						run("Select None");
						selectWindow(background_mask_title); //"background mask"); //02_10_21 SK
						//waitForUser("before Invert 1");
						//run("Invert");
						run("Create Selection");
						//waitForUser("after Invert 1");
						
						selectWindow(list[j]);
						run("Restore Selection");
						
						//loop over all frames...
						mean_bk_arr = newArray(0);
						for (n=1; n<=nSlices; n++) 
						{
	          				setSlice(n);
	          				getRawStatistics(count, mean, min, max, std,histogram);
							mean_bk_arr = Array.concat(mean_bk_arr, mean);
						}
						
						//line
						selectWindow(list[j]);
						run("Select None");
						
						selectWindow(nuclei_line_mask_title); //"nuclei line mask"); //02_10_21 SK
						//waitForUser("before Invert 2");
						//run("Invert");
						run("Create Selection");
						//waitForUser("after Invert 2");
						
						selectWindow(list[j]);
						run("Restore Selection");
						
						//loop over all frames...
						mean_line_arr = newArray(0);
						mean_line_to_bk_arr = newArray(0);
						for (n=1; n<=nSlices; n++) 
						{
	          				setSlice(n);
	          				getRawStatistics(count, mean, min, max, std,histogram);
							mean_line_arr = Array.concat(mean_line_arr, mean);
							mean_line_to_bk_arr = Array.concat(mean_line_to_bk_arr, mean/mean_bk_arr[n-1]);
						}
						
						//plot
						Plot.create("Plot mean line intensity", "frame", "mean intensity");
						Plot.add("line",mean_line_arr);
						Plot.show();
						
						Plot.create("Plot mean backgr intensity", "frame", "mean intensity");
						Plot.add("line",mean_bk_arr);
						Plot.show();
						
						Plot.create("Plot ratio intensity", "frame", "intensity ratio");
						Plot.add("line",mean_line_to_bk_arr);
						Plot.show();
						
						//show nuclei selection for user approval
						selectWindow(list[j]);
						run("Select None");
						selectWindow(background_mask_title); //"background mask"); //02_10_21 SK
						run("Create Selection");
						selectWindow(list[j]);
						run("Restore Selection");
						setSlice(1);
						
						//arrange windows for easy viewing
						selectWindow("Plot mean line intensity");
						getLocationAndSize(x, y, width, height);
						
						setLocation(scr_w*0.2,scr_h*0.1,width/2,height/2);
						getLocationAndSize(x, y, width, height);
						
						selectWindow("Plot mean backgr intensity");
						setLocation(x,y+height,width,height);
						getLocationAndSize(x, y, width, height);
						
						selectWindow("Plot ratio intensity");
						setLocation(x,y+height,width,height);
						
						selectWindow(list[j]);
						setLocation(x+width,y);
						
						Dialog.createNonBlocking("Continue?");
						Dialog.setLocation(x, y+height+height);
					
						Dialog.addMessage("Are the segmentation and plots okay? (Uncheck to remove this data from results.)");
						Dialog.addCheckbox("Good alignment", 1);
						Dialog.addNumber("Fail frame", fail_frame, 0, 3, "");
						Dialog.show();
						retVal=Dialog.getCheckbox();
						fail_frame=Dialog.getNumber();
						
						if(retVal)
						{	
							run("Set Measurements...", "area mean min integrated redirect=None decimal=3");
						
							//background
							//loop over all frames...
							selectWindow(list[j]);
							for (n=1; n<=nSlices; n++) 
							{
	          					setSlice(n);
	          					run("Measure");
	          					
								getRawStatistics(count, mean, min, max, std,histogram);
								
								total_sat=0;
	          					if(max >= PixelMax) //alignment may change pixel values, so possible > 4095
	          					{
	          						for(m=PixelMax; m<=max; m++)
	          						{
	          							total_sat += histogram[m];
	          						}
	          					}
	          					setResult("Frame",nResults-1, n);
	          					setResult("FailFrame", nResults-1, fail_frame);
	          					setResult("Name", nResults-1, list[j]);
	          					setResult("RoiType", nResults-1, "nuclei");
	          					setResult("RatioSaturated", nResults-1, total_sat/count);
	          					updateResults();
	       					}
						
							
						
							//line
							selectWindow(list[j]);
							run("Select None");
						
							selectWindow(nuclei_line_mask_title); //"nuclei line mask"); //02_10_21 SK
							run("Create Selection");
							selectWindow(list[j]);
							run("Restore Selection");

							//loop over all frames...
							for (n=1; n<=nSlices; n++) 
							{
	          					setSlice(n);
	          					run("Measure");
	          					
	          					getRawStatistics(count, mean, min, max, std,histogram);

								total_sat=0;
	          					if(max >= PixelMax) //alignment may change pixel values, so possible > 4095
	          					{
	          						for(m=PixelMax; m<=max; m++)
	          						{
	          							total_sat += histogram[m];
	          						}
	          					}
	          					setResult("Frame",nResults-1, n);
	          					setResult("FailFrame", nResults-1, fail_frame);
	          					setResult("Name", nResults-1, list[j]);
	          					setResult("RoiType", nResults-1, "line");
	          					setResult("RatioSaturated", nResults-1, total_sat/count);
	          					updateResults();
	       					}

							saveAs("Results",all_csv_name);
							close("*");
						}
						else
						{
							File.delete(dir+'/'+roi_selection2);
							close("*");
						}
						
					}
					else 
					{//else: bad alignment, do nothing more
						
						close("*");
					}
				}
			}
		}
		
		run("Clear Results");
		
		//selectWindow("Results");
		//run("Close");
	}
}



function listDirs(dir)
{
     list = getFileList(dir);
	 dirList = newArray();
     for (i=0; i<list.length; i++)
     {
        if (endsWith(list[i], "/"))
        {
        	//print(list[i]);
           	dirList = Array.concat(dirList, list[i]);
        }
     }
     return dirList;
}


