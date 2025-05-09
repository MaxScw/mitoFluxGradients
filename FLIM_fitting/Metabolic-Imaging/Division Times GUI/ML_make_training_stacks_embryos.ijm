
// Paste multiple calls with different paths here to run a batch of segmentations
EmbryosSeg("Z:\\Lab\\Tim\\2016-10-03 Live Birth Acquisitions\\2017-04-04 Batch5\\s1_a1")






function EmbryosSeg(path){
newImage("Untitled", "8-bit black", 512, 512, 1);
run("Trainable Weka Segmentation");
wait(1000)
//selectWindow("Trainable Weka Segmentation v3.2.5");
run("Close All");

// Get dependent paths and load IllProfCal_fad.tif
DailyPath = "\\";
parts1=split(path, "\\");
l = parts1.length;
for(i=0; i<parts1.length-1; i++) { 
DailyPath = DailyPath+parts1[i]+"\\";
}
DailyPath = DailyPath+"DailyFiles\\";

open(DailyPath+"illprofcal_fad.tif");
rename("IllProf");
run("Select All");
getStatistics(Illarea,Illmean);



// Get sub-directories of tiffs from 'MultiD_tiff_convert.m'
sdtpath = path+"\\sorted_sdts\\";
Dsdt = getFileList(sdtpath);
L = Dsdt.length;

//Trained Weka classifier model is stored in 'Daily Files'. 
// Should generally be consistent for all data on a given day.
Classifierpath = DailyPath+"classifier.model";


for(i=0; i<L; i++) { 
	// For Intensity Tiff folders only. Disregard other folders and files.
	if(startsWith(Dsdt[i],"IntTiffs_")){ //&endsWith(Dsdt[i],"\\")
		print(i);
		//waitForUser("fetch, scum");

		
		// Open IllProfCal image
		open(DailyPath+"illprofcal_fad.tif");
		rename("IllProf");
		run("Select All");
		getStatistics(Illarea,Illmean);

		run("Image Sequence...", "open=["+sdtpath+Dsdt[i]+"]");
		rename("FAD");
		makeRectangle(0, 0, 512, 512);
		run("Duplicate...", "duplicate");
		rename("NADH");
		imageCalculator("Divide create 32-bit stack", "NADH","IllProf");
		close("NADH");
		selectWindow("Result of NADH");
		rename("NADH");
		run("Multiply...", "value="+Illmean+" stack");
		
		selectWindow("FAD");
		makeRectangle(512, 0, 512, 512);
		run("Crop");
		imageCalculator("Divide create 32-bit stack", "FAD","IllProf");
		close("FAD");
		selectWindow("Result of FAD");
		rename("FAD");
		run("Multiply...", "value="+Illmean+" stack");

		
		// Multiply NADH and FAD together to get a joint image. Will use this both for 
		// disk convolution (WoW crop) and Weka Segmentation
		imageCalculator("Multiply create 32-bit stack", "NADH","FAD");
		rename("DualMult");
		close("NADH");
		close("FAD");
		//Adjust contrast
		makeOval(106, 106, 300, 300);
		run("Enhance Contrast", "saturated=0.35");
		run("Select None");
		run("8-bit");

		// Crop out WoW dish edge
		run("Duplicate...", "title=CropIm");
		run("Gaussian Blur...", "sigma=10");
		setAutoThreshold("Default dark");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Invert");
		run("Keep Largest Region");
		close("CropIm");
		rename("CropIm");
		run("Fill Holes");
		run("Options...", "iterations=13 count=1 black do=Nothing");
		run("Erode");
		run("Create Selection");
		selectWindow("DualMult");
		run("Restore Selection");
		run("Clear Outside", "stack");
		run("Select None");
		//Save images for training and testing set
		PosL = lengthOf(Dsdt[i]);
		run("Save", "save=["+sdtpath+"ML_"+substring(Dsdt[i], 9, PosL-1)+".tif]");
		close("CropIm");
		
	};
};		
};