/* $id$ */

/*
Adapted after working with EIGER HDF5 to CBF converter
 which was written by Takanori Nakane

imode = 1 - read the number of images in the dataset only
      = 2 - read all the header metadata
      = 3 - read the image data
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <unistd.h>
unsigned int sleep(unsigned int seconds);

extern const H5Z_class2_t H5Z_LZ4;
extern const H5Z_class2_t bshuf_H5Filter;
void register_filters() {
  H5Zregister(&H5Z_LZ4);
  H5Zregister(&bshuf_H5Filter);
}

int fromhdf5_(char *imageFile, unsigned short *image, int *firstimg, int *lastimg, 
	      int *numimages, int cbf_int[16], int cbf_int4[16], 
	      double cbf_double[16],char cbf_char[16][24], double *phibegvals,
	      double *phiendvals, int *imode, int *filetype, int *debug, int *ierr) {
  
  char header[4096] = {};
  int xpixels = -1, ypixels = -1, beamx = -1, beamy = -1, nimages = -1, ntrigger = -1, depth = -1, countrate_cutoff = -1, image_nr_high = -1;
  int from, to, ic,nophivals,mode,ifiletype,idebug;
  int ret;
  int block_start = 1;
  int static pixel_warn = 0;
  int static bit_depth_warn = 0;
  double pixelsize_x = -1, pixelsize_y = -1, wavelength = -1, distance = -1, count_time = -1, frame_time = -1, osc_width = -1, osc_start = -9999, thickness = -1;
  double phirng,sched_start,sched_end;
  char detector_sn[256] = {}, description[256] = {}, version[256] = {}, pixelsize_units[8] = {}, distance_units[8] = {}, dataset[256] = {};
  
  hid_t hdf;
  hid_t entry, group,data,dataspace;
  herr_t status;
  hsize_t ndimensions,dimensions[3];
  hsize_t dims[3];
  char data_name[80] = {};
  char data_set[80] = {};
  int number_per_block = 0;
  int xypixels[2];
  if (strstr(imageFile,".nxs") == NULL && strstr(imageFile,"_master.h5") == NULL) {
    printf("this filename is not a DLS '.nxs' file or a classic HDF5 'master' file - EXITING NOW\n");
    return -2;
  }
  
  /* ------ need to go through this very carefully to remove duplicate code and tidy things up
            but this current code should work okay for classic HDF5 master files and DLS style .NXS 
            files for the moment.

	    Be very careful when changing code - it may be necessary to change for both file types, 
	    or it may be necessary to NOT change the code!
	    -----*/


  /* ---- classic style files with metadata contained in _master.h5 file ---- */
  /*      DLS-style NXS files follow at around ine 541 on 27/03/2018          */

  if (strstr(imageFile,"_master.h5") != NULL) {
  if (idebug != 0 )  printf("hdf5 file is %s \n",imageFile); 

  register_filters();
  from = *firstimg;
  to = *lastimg;
  mode = *imode;
  idebug = *debug;
  ifiletype = *filetype;
  if (idebug != 0 )  printf("mode is %i \n",mode); 
  if (idebug != 0 ) printf("ifiletype is %i \n",ifiletype);

  hdf = H5Fopen(imageFile, H5F_ACC_RDONLY, H5P_DEFAULT); 

  if (idebug != 0 ) printf("hdf is %i \n",hdf);

  if (hdf < 0) {
    fprintf(stdout, "Failed to open file %s\n", imageFile);
    return -1;
  }
  if (idebug != 0 ) printf("ifiletype is %i \n",ifiletype);

  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/nimages", &nimages);
  if (idebug != 0 ) printf("nimages  %d \n",nimages);

  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/ntrigger", &ntrigger);
  if (idebug == 1 ) printf("ntrigger %d \n",ntrigger);

  // Eiger images do not seem to have the field image_nr_high so this causes an HDF5 error (non-fatal)
  //H5LTread_dataset_int(hdf, "/entry/data/data/image_nr_high", &image_nr_high);
  //if (idebug == 1 )  printf("nimages from image_nr_high  %d \n",image_nr_high);

  if ( ntrigger > 1 ) {
      nimages = ntrigger * nimages;
  }

  *numimages = nimages/(1+to-from);
  if (idebug == 1 ) printf("numimages set to  %d \n",*numimages);
  if ( *imode == 3 ) {
    H5Fclose(hdf);
      return 0;
  }


  H5Eset_auto(0, NULL, NULL); // Comment out this line for debugging.
  H5LTread_dataset_string(hdf, "/entry/instrument/detector/description", description);
  H5LTread_dataset_string(hdf, "/entry/instrument/detector/detector_number", detector_sn);
  H5LTread_dataset_string(hdf, "/entry/instrument/detector/detectorSpecific/software_version", version);
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/bit_depth_image", &depth);
  if (depth > 0) {

    // this is okay!

  } else {
    if ( bit_depth_warn == 0 ) {
      bit_depth_warn = 1;
      fprintf(stdout, " WARNING: /entry/instrument/detector/bit_depth_image is not available. We assume 16 bit.\n");
    }
    depth = 16;
  }
  unsigned int error_val = (unsigned int)(((unsigned long long)1 << depth) - 1);
  // Saturation value
  // Firmware >= 1.5
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/saturation_value", &countrate_cutoff);
  if (countrate_cutoff < 0) {
    // Firmware >= 1.4
    H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/countrate_correction_count_cutoff", &countrate_cutoff);
    if (countrate_cutoff > 0) {
      countrate_cutoff++;
    } else {
      // < 1.4
      H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/detectorModule_000/countrate_correction_count_cutoff", &countrate_cutoff);
      if (countrate_cutoff > 0) {
        countrate_cutoff++;
      } else {
        countrate_cutoff = error_val - 1;
      }
    }
  }

  H5LTread_dataset_double(hdf, "/entry/instrument/detector/sensor_thickness", &thickness); // in m 
  if (thickness <= 0) {
    thickness = 450E-6;
  }
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/x_pixels_in_detector", &xpixels);
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/y_pixels_in_detector", &ypixels);
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/beam_center_x", &beamx);
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/beam_center_y", &beamy);
  H5LTread_dataset_double(hdf, "/entry/instrument/detector/count_time", &count_time);
  H5LTread_dataset_double(hdf, "/entry/instrument/detector/frame_time", &frame_time);
  H5LTread_dataset_double(hdf, "/entry/instrument/detector/x_pixel_size", &pixelsize_x);
  H5LTread_dataset_double(hdf, "/entry/instrument/detector/y_pixel_size", &pixelsize_y);

  // Detector distance

  H5LTread_dataset_double(hdf, "/entry/instrument/detector/distance", &distance); // Firmware >= 1.7
  if (idebug == 1 ) printf("distance %12.4f \n",distance);  
  if (distance > 0) {
  } else {

    H5LTread_dataset_double(hdf, "/entry/instrument/detector/detector_distance", &distance); // Firmware< 1.7
    if (idebug == 1 )  printf("distance from detector_distance %12.4f \n",distance); 
  }
  // Test if distance has been defined in mm instead of m, assume detector at least 10mm away
  if (distance > 10) {
      distance = 0.001 * distance;
	 }
  // Wavelength
  H5LTread_dataset_double(hdf, "/entry/sample/beam/incident_wavelength", &wavelength); // Firmware >= 1.7
  if (wavelength <= 0) {

    H5LTread_dataset_double(hdf, "/entry/instrument/beam/wavelength", &wavelength);
    if (wavelength <= 0) {

      H5LTread_dataset_double(hdf, "/entry/instrument/monochromator/wavelength", &wavelength);
      if (wavelength <= 0) {
	H5LTread_dataset_double(hdf, "/entry/instrument/beam/incident_wavelength", &wavelength); // Firmware 1.6
      }
    }
  }
   printf("wavelength %12.5f \n",wavelength); 
  if (wavelength < 0) {
    fprintf(stdout, " WARNING: wavelength was not defined! \"Wavelength\" field in the output is set to -1.\n");
  }
  
  if (idebug == 1 )  printf("nimages just before new read %d \n",nimages);
  /*
    Repeat read of nimages here because for DLS HDF5 master files the value read
    at (or near) line 80, gets reset to zero after the call H5LTread_dataset_string.
    Not at all clear why or how this happens (does not happen with ESRF files)  AGWL 4 Apr 2019
  */
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/nimages", &nimages);
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/ntrigger", &ntrigger);

  if (idebug == 1 ) {
    printf("nimages just after new read %d \n",nimages);
    printf("ntrigger %d \n",ntrigger);
    printf("numimages after new read  %d \n",*numimages);
  }


  if ( ntrigger > 1 ) {
      nimages = ntrigger * nimages;
  }

  
  /*
    check that values exist for phi for ALL images, if not try other locations in metadata
  */
  if (ifiletype == 0) {
  H5LTread_dataset_double(hdf, "/entry/sample/goniometer/phi", phibegvals);
  H5LTread_dataset_double(hdf, "/entry/sample/goniometer/phi_end", phiendvals);
  if (idebug == 1 ) {
     printf("phi start 1st image %8.2f  end %8.2f nimages %d \n",phibegvals[0],phiendvals[0],nimages); 
     printf("phi start 2nd image %8.2f  end %8.2f nimages %d \n",phibegvals[1],phiendvals[1],nimages);
  }
  nophivals = 0;
  /*
    special test if only one image
  */
  if (nimages == 1) {
      if (phibegvals[0] == phiendvals[0]) nophivals =1;
  } else {
      for (ic = 0; ic < nimages; ic++) {
        if (phibegvals[ic] == phiendvals[ic]) nophivals = 1;
      }
  }
  if (nophivals == 0) {
        printf ("rotation angles are stored in goniometer/phi\n");
        *filetype = 1;
  } else {
    nophivals = 0;
    H5LTread_dataset_double(hdf, "/entry/sample/goniometer/omega", phibegvals);
    H5LTread_dataset_double(hdf, "/entry/sample/goniometer/omega_end", phiendvals);
    if (idebug == 1 ) {
      printf("omega start 1st image %8.2f  end %8.2f\n",phibegvals[0],phiendvals[0]); 
      printf("omega start 2nd image %8.2f  end %8.2f\n",phibegvals[1],phiendvals[1]); 
      printf("omega start 3rd image %8.2f  end %8.2f\n",phibegvals[2],phiendvals[2]); 
      printf("omega start 4th image %8.2f  end %8.2f\n",phibegvals[3],phiendvals[3]); 
    }
    /*
      special test if only one image
    */
    if (nimages == 1) {
        if (phibegvals[0] == phiendvals[0]) nophivals =1;
    } else {
        for (ic = 0; ic < nimages; ic++) {
         if (phibegvals[ic] == phiendvals[ic]) nophivals = 1;
        }
    }
    if (nophivals == 0) {
        printf ("rotation angles are stored in goniometer/omega\n");
        *filetype = 2;
    } else {
	// try DLS master file formats
       nophivals = 0;
       H5LTread_dataset_double(hdf, "/entry/sample/transformations/omega", phibegvals);
       H5LTread_dataset_double(hdf, "/entry/sample/transformations/omega_end", phiendvals);
       printf("DLS omega start 1st image %8.2f  end %8.2f\n",phibegvals[0],phiendvals[0]);
       nophivals = 0;
       for (ic = 0; ic < nimages; ic++) {
         if (phibegvals[ic] == phiendvals[ic]) nophivals = 1;
       }
       if (nophivals == 0) {
            printf("rotation angles are in transformations/omega \n");
            *filetype = 3;
       }
    }
    }
  H5LTread_dataset_double(hdf, "/entry/sample/goniometer/omega_range_average", &phirng);
    if (phirng < 0) {
      phirng = 0;
    }
    if (phirng == (phiendvals[nimages-1] - phibegvals[0]) / nimages) {
	if (idebug == 1 ) printf("calculated and stored phi range match: %8.2f\n",phirng);
    } else {
	if (idebug == 1 )  printf("stored phi range is %8.2f, calculated is %8.2f AND they do not match!\n",phirng,(phiendvals[nimages-1] - phibegvals[0]) / nimages);
    }
    

  H5LTread_dataset_double(hdf, "/entry/sample/goniometer/omega_range_average", &osc_width);
  if (osc_width > 0) {
      if (idebug == 1 )  fprintf(stdout, " /entry/sample/goniometer/omega_range_average = %f (deg)\n", osc_width);
  } else {
      if (idebug == 1 )  fprintf(stdout, " WARNING: oscillation width was not defined. \"Start_angle\" field in the output is set to 0!\n");
    osc_width = 0;
  }
  }

  int bufsize = (nimages < 100000) ? 100000 : nimages;
  double *angles = (double*)malloc(bufsize * sizeof(double));
  /* Next print previously include to avoid loosing values of "from,to"
     in DLS HDF5 files when building on CentOS */
  printf("After malloc call, from %d to %d \n",from,to);
  // I don't know why but nimages can be too small ...
  if (angles == NULL) {
    fprintf(stdout, "failed to allocate buffer for omega.\n");
    return -1;
  }
  angles[0] = -9999;
  H5LTread_dataset_double(hdf, "/entry/sample/goniometer/omega", angles);

  entry = H5Gopen2(hdf, "/entry", H5P_DEFAULT);
  if (entry < 0) {
    fprintf(stdout, "/entry does not exist!\n");
    return -1;
  }

  strcpy(&cbf_char[5][0],description);
  strcat(&cbf_char[5][0],"\0");
  strcpy(&cbf_char[7][0],detector_sn);
  strcat(&cbf_char[7][0],"\0");
  strcpy(&cbf_char[8][0],version);
  strcat(&cbf_char[8][0],"\0");

  /* ---- output section & update cbf_char, cbf_int, etc ----- */
  depth = 16;

  error_val = (unsigned int)(((unsigned long long)1 << depth) - 1);
  cbf_int4[0] = countrate_cutoff;
  cbf_int[2] = xpixels;
  cbf_int[3] = ypixels;
  cbf_int[4] = beamx;
  cbf_int[5] = beamy;
  cbf_int[6] = 0;
  cbf_int[7] = 32767;

  /* Need to re-assign "from,to" here otherwise CentOS build does not
     read DLS files correctly */

  from = *firstimg;
  to = *lastimg;
  if (idebug == 1 ) {
    printf("from is %i \n",from); 
    printf("to is %i \n",to); 
    printf("numimages is now  %d \n",*numimages);
  }

  cbf_double[0] = wavelength;
  cbf_double[1] = distance * 1000.0;
  cbf_double[2] = pixelsize_x * 1000.0;
  cbf_double[3] = pixelsize_y * 1000.0;
  cbf_double[5] = phibegvals[from-1];
  cbf_double[6] = phiendvals[to-1];

  if (idebug == 1 ) {
    printf("cbf_double[5] %12.4f \n",cbf_double[5]);
    printf("cbf_double[6] %12.4f \n",cbf_double[6]);
  }

  cbf_double[7] = cbf_double[6] - cbf_double[5];
  cbf_double[11] = thickness * 1000.0;
  cbf_double[12] = 0.0;

  printf("finished reading header, mallocing buffer with size %d x %d\n",xpixels,ypixels);
  fflush (stderr);
  fflush (stdout);

  signed int *multibuf = (signed int*)malloc(sizeof(signed int) * xpixels * ypixels);
  unsigned int *buf = (unsigned int*)malloc(sizeof(unsigned int) * xpixels * ypixels);
  signed int *buf_signed = (signed int*)malloc(sizeof(signed int) * xpixels * ypixels);
  signed int *pixel_mask = (signed int*)malloc(sizeof(signed int) * xpixels * ypixels);
  for (ic = 0; ic < xpixels * ypixels; ic++) {
    multibuf[ic] = 0;
  }
  if (multibuf == NULL || buf == NULL || buf_signed == NULL || pixel_mask == NULL) {
    fprintf(stdout, "Failed to allocate image buffer.\n");
    return -1;
  }

  // Check if /entry/data present
  group = H5Gopen2(hdf, "/entry/data", H5P_DEFAULT); 
  if (group < 0) {
    group = entry; // leak!
  } 
  if (H5LTfind_dataset(group, "data_000000")) {
      block_start = 0;
  } else {
      if (idebug == 1 )  fprintf(stdout, "This dataset starts from data_000001.\n");
      block_start = 1;
  }

  if (idebug == 1) printf("block_start is:  %d \n",block_start);
  fflush (stdout);

  snprintf(data_name, 20, "data_%06d", block_start);
  

  pixel_mask[0] = -9999;
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/pixel_mask", pixel_mask);

  if (pixel_mask[0] == -9999 && pixel_warn == 0) {
    pixel_warn = 1;
    fprintf(stdout, "WARNING: failed to read the pixel mask from /entry/instrument/detector/detectorSpecific/pixel_mask.\n");
    fprintf(stdout, " Thus, we mask pixels whose intensity is %u (= (2 ^ bit_depth_image) - 1) by converting them to -1. \n", error_val);
    fprintf(stdout, " However, this might mask overloaded (saturated) pixels as well.\n");
  }

  fflush (stdout);
 
  
  // Check if /entry/data present
  group = H5Gopen2(entry, "data", H5P_DEFAULT);  
  if (group < 0) {
    group = entry; // leak!
  }

  number_per_block = 0;
  
  // Open the first data block to get the number of frames in a block
  if (idebug == 1) printf("block_start is now:  %d \n",block_start);
  snprintf(data_name, 20, "data_%06d", block_start);
  if (idebug == 1) printf("data_name is:  %s \n",data_name);
   
  data = H5Dopen2(group, data_name, H5P_DEFAULT);
  if (idebug == 1) printf("data is:  %d \n",data);
  if (data < 0) {
    fprintf(stdout, "ierr=10, failed to open /entry/%s\n", data_name);
    printf(stdout, "ierr=10, failed to open /entry/%s\n", data_name);
    if (idebug == 1) printf("numimages is now  %d \n",*numimages);
    *ierr = 10;
    fflush (stdout);
    return 0;
  }
  dataspace = H5Dget_space(data);
  if (H5Sget_simple_extent_ndims(dataspace) != 3) {
    fprintf(stdout, "Dimension of /entry/%s is not 3!\n", data_name);
    return -1;    
  }
  H5Sget_simple_extent_dims(dataspace, dims, NULL);
  number_per_block = dims[0];
  //  fprintf(stdout, "The number of images per data block is %d.\n", number_per_block);

  H5Sclose(dataspace);
  H5Dclose(data);

  int frame;

  for (frame = from; frame <= to; frame++) {

    if (angles[0] != -9999) {
      osc_start = angles[frame - 1];
    } else {
      osc_start = osc_width * frame; // old firmware
    }
    
    if (frame > nimages) {
      // Due to a firmware bug, nimages can be smaller than the actual value.
      // So we don't exit here
    }
    
    int block_number = block_start + (frame - 1) / number_per_block;
    int frame_in_block = (frame - 1) % number_per_block;
    snprintf(data_name, 20, "data_%06d", block_number); 
    data = H5Dopen2(group, data_name, H5P_DEFAULT);
    dataspace = H5Dget_space(data);
    if (data < 0) {
      fprintf(stdout, "failed to open /entry/%s\n", data_name);
      *ierr = 10;
      return 0;
    }
    if (H5Sget_simple_extent_ndims(dataspace) != 3) {
      fprintf(stdout, "Dimension of /entry/%s is not 3!\n", data_name);
      return -1;    
    }
    
    // Get the frame
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    hsize_t offset_in[3] = {frame_in_block, 0, 0};
    hsize_t offset_out[3] = {0, 0, 0};
    hsize_t count[3] = {1, ypixels, xpixels};
    hid_t memspace = H5Screate_simple(3, dims, NULL); 
    if (memspace < 0) {
      fprintf(stdout, "failed to create memspace\n");
      return -1;
    }
    
    ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset_in, NULL, 
			      count, NULL);
    if (ret < 0) {
      fprintf(stdout, "select_hyperslab for file failed\n");
      *ierr = 20;
      return 0;
    }
    ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, 
			      count, NULL);
    if (ret < 0) {
      fprintf(stdout, "select_hyperslab for memory failed\n");
      *ierr = 20;
      return 0;
    }
    ret = H5Dread(data, H5T_NATIVE_UINT, memspace, dataspace, H5P_DEFAULT, buf);
    if (ret < 0) {
      fprintf(stdout, "H5Dread for image failed. Wrong frame number?\n");
      *ierr = 20;
      return 0;
    }
    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Dclose(data);
    for (ic = 0; ic < xpixels * ypixels; ic++) {
      if ( buf[ic] == -1 || multibuf[ic] == -1 ) {
	multibuf[ic] = -1;
      } else if ((pixel_mask[0] != -9999 && pixel_mask[ic] == 1) ||  // the pixel mask is available
		 (pixel_mask[0] == -9999 && buf[ic] == error_val)) {  // not available or
	multibuf[ic] = -1 ;
      } else if ((pixel_mask[0] != -9999 && pixel_mask[ic] > 1) ||   // the pixel mask is 2, 4, 8, 16 or
		 (multibuf[ic] == -2 )) {                            // already set 
	multibuf[ic] = -2 ;
      } else {
	multibuf[ic] = multibuf[ic] + buf[ic];
	cbf_double[12] += (double) buf[ic];
	
	// set up max & min counts per pixel on the summed image
	cbf_int[6] = cbf_int[6] > multibuf[ic] ? cbf_int[6] : multibuf[ic];
	cbf_int[7] = cbf_int[7] < multibuf[ic] ? cbf_int[7] : multibuf[ic];
      }
    }

    
    ctimer_(&sched_end);

    sched_start = sched_end;

    for (ic = 0; ic < xpixels * ypixels; ic++) {
      /*
	superpacking from Andrew- transfer to Mosflm IMAGE array
      */
      
      int excount = 0;
      int hicount = 0;
      int midcount = 0;
      int lowcount = 0;
      int undercount = 0;
      //    if (multibuf[ic] <= 1000)lowcount++;
      //    if (multibuf[ic] > 1000 && multibuf[ic] <= 8000 )midcount++;
      //    if (multibuf[ic] > 8000 && multibuf[ic] <= 16000)hicount++;
      //    if (multibuf[ic] > 16000)excount++;
      if (multibuf[ic] < 0) {
	undercount++;
      }
      if (multibuf[ic] <= 32767) {
	image[ic] = (short) multibuf[ic];
	lowcount++;
      }
      else if (multibuf[ic] < 131072) {
	image[ic] = (short)  -((multibuf[ic]+8)/16) + 2045;
	midcount++;
      }
      else if (multibuf[ic] < 3538000) {
	image[ic] = (short) -((multibuf[ic]+64)/128) - 5124;
	hicount++;
      }
      else {
	image[ic] = (short) -32766;
	excount++;
      }
      lowcount -= undercount;
      image[ic] = (short) multibuf[ic];
    }
  }
    /* --------- close group & file -------- */
    H5Gclose(group);
    H5Fclose(hdf);
      


    ////// free memory -----
    free(multibuf);
    free(buf);
    free(buf_signed);
    free(angles);
    free(pixel_mask);

      
      return 0;
  
}

  /* ---- DLS style files with metadata contained in .NXS file ---- */

  if (strstr(imageFile,".nxs") != NULL) {

  register_filters();
  
  from = *firstimg;
  to = *lastimg;
  
  hdf = H5Fopen(imageFile, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (hdf < 0) {
    fprintf(stdout, "Failed to open file %s\n", imageFile);
    return -1;
  }
  
    status = H5Gget_objinfo (hdf, "/entry/instrument/eiger", 0, NULL);
    if (status == 0 ) {
      entry = H5Gopen2(hdf, "/entry/instrument/eiger", H5P_DEFAULT);
      if (entry < 0) {
	fprintf(stdout, "/entry does not exist!\n");
	return -1;
      }
    }


    snprintf(data_name, 20, "/entry/eiger/omega"); 
    group = H5Dopen2(hdf, data_name, H5P_DEFAULT);


    dataspace = H5Dget_space(group);

    H5Sget_simple_extent_dims(dataspace,dimensions,NULL);
    nimages = (int) dimensions[0];
    *numimages = nimages/(1+to-from);
    
    H5Dclose(group);

    // H5Eset_auto(NULL, NULL); // Comment out this line for debugging.
    
    /*-------------------------------
    need to get image dimensions somehow....

    it's stored in dataspace /entry/instrument/eiger/module/data_size
    __________________________________*/
    //    xpixels = 2070;
    //    ypixels = 2167;
    // printf("images are %d pixels by %d pixels\n",xpixels,ypixels);	
    H5LTread_dataset_int(hdf, "/entry/instrument/eiger/module/data_size",xypixels);
    xpixels = xypixels[0];
    ypixels = xypixels[1];
    /* -----------------------------*/
    /* detector details */
    H5LTread_dataset_string(hdf, "/entry/instrument/eiger/description", description);

    /* DLS files don't have the serial number or software version */
    strcpy(detector_sn,"unknown");
    strcpy(version,"unknown");
    H5Eset_auto(0, NULL, NULL); // Comment out this line for debugging.
    // Saturation value
    // Firmware >= 1.5
    H5LTread_dataset_int(hdf, "/entry/instrument/eiger/saturation_value", &countrate_cutoff);
    H5LTread_dataset_double(hdf, "/entry/instrument/eiger/sensor_thickness", &thickness); // in m
    H5LTread_dataset_int(hdf, "/entry/instrument/eiger/beam_centre_x", &beamx);
    H5LTread_dataset_int(hdf, "/entry/instrument/eiger/beam_centre_y", &beamy);
    H5LTread_dataset_double(hdf, "/entry/instrument/eiger/count_time", &count_time);
    H5LTread_dataset_double(hdf, "/entry/instrument/eiger/x_pixel_size", &pixelsize_x);
    H5LTread_dataset_double(hdf, "/entry/instrument/eiger/y_pixel_size", &pixelsize_y);  
    H5LTread_dataset_double(hdf, "/entry/instrument/detector_z/det_z", &distance);
    H5LTread_dataset_double(hdf, "/entry/sample/beam/incident_wavelength", &wavelength);
    // DLS NXS files only have phi start....
    H5LTread_dataset_double(hdf, "/entry/eiger/omega", phibegvals);
    phirng = phibegvals[1] - phibegvals[0];
    for (ic = 0; ic < nimages-1; ic++) {
      phiendvals[ic] = phibegvals[ic] + phirng;
    }
    
    snprintf(data_name, 20, "data_0"); 

  int bufsize = (nimages < 100000) ? 100000 : nimages;
  double *angles = (double*)malloc(bufsize * sizeof(double));
  // I don't know why but nimages can be too small ...
  if (angles == NULL) {
    fprintf(stdout, "failed to allocate buffer for omega.\n");
    return -1;
  }
  angles[0] = -9999;
  if (strstr(imageFile,".nxs") != NULL) {
    H5LTread_dataset_double(hdf, "/entry/eiger/omega", angles);
  }
  if (strstr(imageFile,"_master.h5") != NULL) {
    H5LTread_dataset_double(hdf, "/entry/sample/goniometer/omega", angles);
  }

 
  strcpy(&cbf_char[5][0],description);
  strcat(&cbf_char[5][0],"\0");
  strcpy(&cbf_char[7][0],detector_sn);
  strcat(&cbf_char[7][0],"\0");
  strcpy(&cbf_char[8][0],version);
  strcat(&cbf_char[8][0],"\0");

  /* ---- output section & update cbf_char, cbf_int, etc ----- */
  depth = 16;

  unsigned int error_val = (unsigned int)(((unsigned long long)1 << depth) - 1);

  if (countrate_cutoff > 0) {
    cbf_int4[0] = countrate_cutoff;
  }
  if (thickness < 0) {
    thickness = 450E-6;
    cbf_double[11] = thickness;
  }
  cbf_int[2] = xpixels;
  cbf_int[3] = ypixels;

  cbf_int[4] = beamx;
  cbf_int[5] = beamy;


  if (strcmp(pixelsize_units,"m")) {
    cbf_double[2] = pixelsize_x;
  } else if (strcmp(pixelsize_units,"mm" )) {
    cbf_double[2] = pixelsize_x * 1000.0;
  }
  if (strcmp(pixelsize_units,"m")) {
    cbf_double[3] = pixelsize_y;
  } else if (strcmp(pixelsize_units,"mm" )) {
    cbf_double[3] = pixelsize_y * 1000.0;
  }

  if (distance > 0) {
    cbf_double[1] = distance;
  } 
  if (wavelength > 0) {
    cbf_double[0] = wavelength;
  } 
  if (phirng > 0){
    cbf_double[5] = phibegvals[from-1];
    cbf_double[6] = phiendvals[to-1];
    cbf_double[7] = cbf_double[6] - cbf_double[5];
    phirng = cbf_double[7];
    osc_width = phirng;
  }

  signed   int *multibuf = (signed int*)malloc(sizeof(signed int) * xpixels * ypixels);
  unsigned int *buf = (unsigned int*)malloc(sizeof(unsigned int) * xpixels * ypixels);
  signed   int *buf_signed = (signed int*)malloc(sizeof(signed int) * xpixels * ypixels);
  signed   int *pixel_mask = (signed int*)malloc(sizeof(signed int) * xpixels * ypixels);
  if (multibuf == NULL || buf == NULL || buf_signed == NULL || pixel_mask == NULL) {
    fprintf(stdout, "Failed to allocate image buffer.\n");
    return -1;
  }
  for (ic = 0; ic < xpixels * ypixels; ic++) {
    multibuf[ic] = 0;
  }

  // Check if /entry/data present


  group = H5Gopen(entry, "/entry/instrument/eiger", H5P_DEFAULT);  
  if (group < 0) {
    group = entry; // leak!
  }
  if (H5LTfind_dataset(group, "data_0")) {
    //        fprintf(stdout,"Looks like an Eiger at DLS\n");
  }
  pixel_mask[0] = -9999;
  if ((int) H5LTfind_dataset(group, "pixel_mask") > 0 ) {
	   H5LTread_dataset_int(group, "pixel_mask", pixel_mask);
  } 
  /*-----
  //
  // don't understand why I can't read the pixel mask even though it exists..
  //
  if (pixel_mask[0] == -9999) {
    fprintf(stdout, "WARNING: failed to read the pixel mask from /entry/instrument/eiger/pixel_mask.\n");
    fprintf(stdout, " Thus, we mask pixels whose intensity is %u (= (2 ^ bit_depth_image) - 1) by converting them to -1. \n", error_val);
    fprintf(stdout, " However, this might mask overloaded (saturated) pixels as well.\n");
  } 
  ----*/


  // Open the first data block to get the number of frames in a block
  data = H5Dopen2(group, data_name, H5P_DEFAULT);
  dataspace = H5Dget_space(data);
  if (data < 0) {
    fprintf(stdout, "failed to open /entry/%s\n", data_name);
    return -1;
  }
  if (H5Sget_simple_extent_ndims(dataspace) != 3) {
    fprintf(stdout, "Dimension of /entry/%s is not 3!\n", data_name);
    return -1;    
  }
  //  hsize_t dims[3];
  H5Sget_simple_extent_dims(dataspace, dims, NULL);
  number_per_block = dims[0];
  //  fprintf(stdout, "The number of images per data block is %d.\n", number_per_block);

  H5Sclose(dataspace);
  H5Dclose(data);

  cbf_double[12] = 0.0;
  int frame;
  for (frame = from; frame <= to; frame++) {
    if (angles[0] != -9999) {
      osc_start = angles[frame - 1];
    } else {
      osc_start = osc_width * frame; // old firmware
    }

    if (frame > nimages) {
      fprintf(stdout, "WARNING: invalid frame number specified. %d is bigger than nimages (%d)\n", frame, nimages);
      // Due to a firmware bug, nimages can be smaller than the actual value.
      // So we don't exit here
    }

    int block_number = block_start + (frame - 1) / number_per_block;
    int frame_in_block = (frame - 1) % number_per_block;
    

    data = H5Dopen2(group, data_name, H5P_DEFAULT);

    dataspace = H5Dget_space(data);
    if (data < 0) {
      fprintf(stdout, "failed to open /entry/%s\n", data_name);
      return -1;
    }
    if (H5Sget_simple_extent_ndims(dataspace) != 3) {
      fprintf(stdout, "Dimension of /entry/%s is not 3!\n", data_name);
      return -1;    
    }

    // Get the frame
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    hsize_t offset_in[3] = {frame_in_block, 0, 0};
    hsize_t offset_out[3] = {0, 0, 0};
    hsize_t count[3] = {1, ypixels, xpixels};
    hid_t memspace = H5Screate_simple(3, dims, NULL); 
    if (memspace < 0) {
      fprintf(stdout, "failed to create memspace\n");
      return -1;
    }

    ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset_in, NULL, 
			      count, NULL);
    if (ret < 0) {
      fprintf(stdout, "select_hyperslab for file failed\n");
      return -1;
    }
    ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, 
			      count, NULL);
    if (ret < 0) {
      fprintf(stdout, "select_hyperslab for memory failed\n");
      return -1;
    }
    ret = H5Dread(data, H5T_NATIVE_UINT, memspace, dataspace, H5P_DEFAULT, buf);
    if (ret < 0) {
      fprintf(stdout, "H5Dread for image failed. Wrong frame number?\n");
      return -1;
    }
    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Dclose(data);



    /////////////////////////////////////////////////////////////////
      cbf_double[0] = wavelength;
      cbf_int[6] = 0;
      cbf_int[7] = 32767;
      for (ic = 0; ic < xpixels * ypixels; ic++) {
	if ( buf[ic] == -1 || multibuf[ic] == -1 ) {
	  multibuf[ic] = -1;
	} else if ((pixel_mask[0] != -9999 && pixel_mask[ic] == 1) ||  // the pixel mask is available
		   (pixel_mask[0] == -9999 && buf[ic] == error_val)) {  // not available or
	  multibuf[ic] = -1 ;
	} else if ((pixel_mask[0] != -9999 && pixel_mask[ic] > 1) ||   // the pixel mask is 2, 4, 8, 16 or
		   (multibuf[ic] == -2 )) {                            // already set 
	  multibuf[ic] = -2 ;
	} else {
	  multibuf[ic] = multibuf[ic] + buf[ic];
	  cbf_double[12] += (double) buf[ic];
	  
	  // set up max & min counts per pixel on the summed image
	  cbf_int[6] = cbf_int[6] > multibuf[ic] ? cbf_int[6] : multibuf[ic];
	  cbf_int[7] = cbf_int[7] < multibuf[ic] ? cbf_int[7] : multibuf[ic];
	}
      }
      
      ctimer_(&sched_end);
      //printf("wrote individual frame to multibuf in %f seconds\n",sched_end-sched_start);
      sched_start = sched_end;
      for (ic = 0; ic < xpixels * ypixels; ic++) {
	/*
	  superpacking from Andrew- transfer to Mosflm IMAGE array
	*/
	
	int excount = 0;
	int hicount = 0;
	int midcount = 0;
	int lowcount = 0;
	int undercount = 0;
	//    if (multibuf[ic] <= 1000)lowcount++;
	//    if (multibuf[ic] > 1000 && multibuf[ic] <= 8000 )midcount++;
	//    if (multibuf[ic] > 8000 && multibuf[ic] <= 16000)hicount++;
	//    if (multibuf[ic] > 16000)excount++;
	if (multibuf[ic] < 0) {
	  undercount++;
	}
	if (multibuf[ic] <= 32767) {
	  image[ic] = (short) multibuf[ic];
	  lowcount++;
	}
	else if (multibuf[ic] < 131072) {
	  image[ic] = (short)  -((multibuf[ic]+8)/16) + 2045;
	  midcount++;
	}
	else if (multibuf[ic] < 3538000) {
	  image[ic] = (short) -((multibuf[ic]+64)/128) - 5124;
	  hicount++;
	}
	else {
	  image[ic] = (short) -32766;
	  excount++;
	}
	lowcount -= undercount;
	image[ic] = (short) multibuf[ic];
      }
  }
      


  /* --------- close group & file -------- */
  H5Gclose(group);
  H5Fclose(hdf);

  ////// free memory -----
  free(multibuf);
  free(buf);
  free(buf_signed);
  free(pixel_mask);
  free(angles);
  
  return 0;
}
  return 0;
}
