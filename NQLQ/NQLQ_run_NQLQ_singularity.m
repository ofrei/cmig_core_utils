function NQLQ_run_NQLQ_singularity(dirname_t1,dirname_flair,dirname_out,sex,age,forceflag,metadata,nq4sif)

if ~exist('forceflag','var') || isempty(forceflag)
  forceflag = false;
  %  forceflag = true;
end
if ~exist('nq4sif','var') || isempty(nq4sif)
  nq4sif='~dale/NQ_Singularity/nq4.sif';
end

starttime0 = now;

corrflag = true;

try
  load('sample_dicom_hdrs.mat');
  fname_complete = sprintf('%s/processing_complete.mat',dirname_out);
  starttime = now();
  if forceflag & exist(dirname_out,'dir')
    cmd = sprintf('rm -rf %s',dirname_out);
    disp(cmd)
    [s o e] = jsystem(sprintf('rm -rf %s',dirname_out)); disp(s); disp(o); disp(e);
    if exist(dirname_out,'dir')
      fprintf('directory %s already exists\n',dirname_out);
    end
  end
  if ~exist(dirname_out,'dir')
    cmd = sprintf('mkdir -p %s',dirname_out);
    disp(cmd)
    [s o e] = jsystem(cmd); disp(s); disp(o); disp(e);
  end
  if corrflag % Perform distortion correction?
    outputdir = sprintf('%s/t1_corr',dirname_out);
    fname_t1_info = sprintf('%s/t1_info.mat',dirname_out);
    if ~exist(fname_t1_info,'file') || forceflag
      if isfile(dirname_t1) % Check if is file, not directory -- assume it's NIFTI
        vol = niftiread_ctx(dirname_t1);
        gradwarpinfo = metadata.gradwarpinfo;
      else
        flist = dir(sprintf('%s/**/*.dcm',dirname_t1)); hdrs = cell([1 length(flist)]); instancevec = NaN([1 length(flist)]); normvec = false([1 length(flist)]); sernumvec = NaN([1 length(flist)]); siuidvec = cell([1 length(flist)]);
        for fi = 1:length(flist)
          fnames{fi} = sprintf('%s/%s',flist(fi).folder,flist(fi).name);
          hdrs{fi} = dicominfo(fnames{fi});
          siuidvec{fi} = hdrs{fi}.SeriesInstanceUID;
          instancevec(fi) = hdrs{fi}.InstanceNumber;
          normvec(fi) = ~isempty(strfind(hdrs{fi}.ImageType,'NORM'));
          sernumvec(fi) = hdrs{fi}.SeriesNumber;
        end
        if length(unique(instancevec))~=length(flist)
          fprintf('directory %s contains more than one series -- choosing normalized one\n',dirname_t1); % Need version that works with GE and Philips
          flist = flist(normvec); hdrs = hdrs(normvec); fnames = fnames(normvec); instancevec = instancevec(normvec);
        end
        [sv si] = sort(instancevec);
        fnames = fnames(si); hdrs = hdrs(si);
        gradwarpinfo = mmil_get_gradwarpinfo(hdrs{1}); [vol_tmp M] = mmil_read_dicom_vol(fnames); vol = ctx_mgh2ctx(vol_tmp,M);
        if ~exist('metadata','var') || isempty(metadata)
          if regexpi(hdrs{1}.Manufacturer,'siemens') 
            metadata.Manufacturer = 'siemens',
          elseif regexpi(hdrs{1}.Manufacturer,'philips') 
            metadata.Manufacturer = 'philips',
          else
            metadata.Manufacturer = 'ge medical';
          end
        end
      end
      switch NQLQ_get_manufacturer(metadata.Manufacturer)
        case 'philips', sample_dicom_hdr = sample_dicom_hdrs_philips.hdr_t1;
        case 'siemens', sample_dicom_hdr = sample_dicom_hdrs_siemens.hdr_t1;
        case 'ge medical', sample_dicom_hdr = sample_dicom_hdrs_ge.hdr_t1;
      end

      nslices = size(vol.imgs,3);
      hdr = sample_dicom_hdr;
      hdr.PatientName = 'Dummy';
      hdr.PatientID = 'Dummy';
      hdrs = repmat({hdr},[1 nslices]);
      Mvxl2lph = vol.Mvxl2lph;
      dcr = Mvxl2lph(1:3,1)/norm(Mvxl2lph(:,1));
        dcc = Mvxl2lph(1:3,2)/norm(Mvxl2lph(:,2));
      dcs = Mvxl2lph(1:3,3)/norm(Mvxl2lph(:,3));
      st = norm(Mvxl2lph(:,3));
      PixelSpacing = colvec(sqrt(sum(Mvxl2lph(:,[1 2]).^2)));
      ImageOrientationPatient = [Mvxl2lph(1:3,1)/norm(Mvxl2lph(1:3,1)); Mvxl2lph(1:3,2)/norm(Mvxl2lph(1:3,2))];
      ImagePositionPatient = Mvxl2lph(1:3,:)*[1 1 1 1]';
      for fi = 1:nslices
% Should also force PatientAge and PatentSex
        hdrs{fi}.SliceThickness = st;
        hdrs{fi}.SpacingBetweenSlices = st;
        hdrs{fi}.PixelSpacing = PixelSpacing;
        hdrs{fi}.ImageOrientationPatient = ImageOrientationPatient;
        hdrs{fi}.ImagePositionPatient = ImagePositionPatient + (fi-1)*dcs*st;
        hdrs{fi}.InstanceNumber = fi;
        hdrs{fi}.SliceLocation = (fi-1)*st;
      end

      if isempty(gradwarpinfo) || ~isfield(gradwarpinfo,'gwtype') % Should also check if correction already performed (gradwarpinfo.unwarpflag)? -- coordinate with Don & Chris C!
        vol_corr = vol;
      else
        vol_corr = ctx_unwarp_grad(vol,gradwarpinfo.gwtype,gradwarpinfo.unwarpflag,gradwarpinfo.isoctrflag);
      end
      [vol_tmp M_tmp] = QD_ctx2mgh(vol_corr); M_tmp = M_RAS_TO_LPH*M_tmp; dims_tmp = size(vol_tmp); ctr = M_tmp(1:3,:)*([(dims_tmp+1)/2 1]');

      if ~exist(outputdir,'dir'), mkdir(outputdir); end

      for fi = 1:length(hdrs)
        if exist('age') && ~isempty(age), hdr.PatientAge = sprintf('%03dY',age); hdr.PatientSex = sex; end
        if isfield(hdrs{fi},'ImageType')
          hdrs{fi}.ImageType = strrep(hdrs{fi}.ImageType,'ND','DIS3D');
        end
        if isfield(hdrs{fi},'ManufacturerModelName')
          hdrs{fi}.ManufacturerModelName = 'Dummy';
        end
        if isfield(hdrs{fi},'ManufacturersModelName')
          hdrs{fi}.ManufacturersModelName = 'Dummy';
        end
        if isfield(hdrs{fi},'AccessionNumber') && length(hdrs{fi}.AccessionNumber)>16
          hdrs{fi}.AccessionNumber = hdrs{fi}.AccessionNumber(1:16);
        end
      end

      hdrs_bak = hdrs;

      % Petend that T1 volume acquired at isocenter, adjust FLAIR accordingly
      if isfield(hdrs{1},'PerFrameFunctionalGroupsSequence') && ~isempty(hdrs{1}.PerFrameFunctionalGroupsSequence.Item_1) % Enhanced dicom?
        if ~isfield(hdrs{fi},'PulseSequenceName'), hdrs{fi}.PulseSequenceName = 't1_mprage'; end
        for slicenum = 1:nslices
          itemname = sprintf('Item_%d',slicenum);
          structname = sprintf('hdrs{fi}.PerFrameFunctionalGroupsSequence.%s.PlanePositionSequence.Item_1',itemname);
          try
            pos = eval(sprintf('%s.ImagePositionPatient',structname))+ctr;
            cmd = sprintf('%s = setfield(%s,''ImagePositionPatient'',pos);',structname,structname); disp(cmd); eval(cmd);
          catch ME
            fprintf(1,'slicenum=%d/%d %s\n',slicenum,nslices,structname);
          end 
        end
      else
        for fi = 1:length(hdrs)
          hdrs{fi}.ImagePositionPatient(:) = hdrs_bak{fi}.ImagePositionPatient(:)-ctr; % This should be hdrs{1}.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1 for enhanced dicom
        end
      end

      dicomwriteVol(vol_corr, outputdir, hdrs); % Should benchmark time used to write DICOM files -- ~1s per image currently

      save(fname_t1_info,'ctr','hdrs');
      fprintf(1,'file %s written\n',fname_t1_info);
    else
      fprintf(1,'file %s already exists\n',fname_t1_info);
      load(fname_t1_info,'ctr');
    end
    dirname_t1_out = outputdir;

    if 0
      vol_t1_dcm = dicomreadVolume_amd(dirname_t1_out); showVol(vol_t1_dcm);
    end

%    vol_corr_t1 = dicomreadVolume_amd(dirname_t1_out); % Doesn't seem to work with enhanced dicom -- update dicomreadVolume_amd to use mmil_read_dicom_vol
    if ~isempty(dirname_flair)
      outputdir = sprintf('%s/flair_corr',dirname_out);
      fname_flair_info = sprintf('%s/flair_info.mat',dirname_out);
      if ~exist(fname_flair_info,'file') || forceflag
        if isfile(dirname_flair) % Check if is file, not directory -- assume it's NIFTI
          vol = niftiread_ctx(dirname_flair);
          gradwarpinfo = metadata.gradwarpinfo;
        else
          flist = dir(sprintf('%s/**/*.dcm',dirname_flair)); hdrs = cell([1 length(flist)]); instancevec = NaN([1 length(flist)]); normvec = false([1 length(flist)]);
          for fi = 1:length(flist)
            fnames{fi} = sprintf('%s/%s',flist(fi).folder,flist(fi).name);
            hdrs{fi} = dicominfo(fnames{fi});
            instancevec(fi) = hdrs{fi}.InstanceNumber;
            normvec(fi) = ~isempty(strfind(hdrs{fi}.ImageType,'NORM'));
          end
          if length(unique(instancevec))~=length(flist)
            fprintf('directory %s contains more than one series -- choosing normalized one\n',dirname_t1); % Need version that works with GE and Philips
            flist = flist(normvec); hdrs = hdrs(normvec); fnames = fnames(normvec); instancevec = instancevec(normvec);
          end
          [sv si] = sort(instancevec);
          fnames = fnames(si); hdrs = hdrs(si);
          gradwarpinfo = mmil_get_gradwarpinfo(hdrs{1}); [vol_tmp M] = mmil_read_dicom_vol(fnames); vol = ctx_mgh2ctx(vol_tmp,M);
        end
        if norm(vol.Mvxl2lph(:,3))>1.5, fieldname = 'hdr_flair_2d'; else fieldname = 'hdr_flair_3d'; end
        switch NQLQ_get_manufacturer(metadata.Manufacturer)
          case 'philips', sample_dicom_hdr = getfield(sample_dicom_hdrs_philips,fieldname);
          case 'siemens', sample_dicom_hdr = getfield(sample_dicom_hdrs_siemens,fieldname);;
          case 'ge medical', sample_dicom_hdr = getfield(sample_dicom_hdrs_ge,fieldname);;
        end

        nslices = size(vol.imgs,3);
        hdr = sample_dicom_hdr;
        hdr.PatientName = 'Dummy';
        hdr.PatientID = 'Dummy';
        hdrs = repmat({hdr},[1 nslices]);
        Mvxl2lph = vol.Mvxl2lph;
        dcr = Mvxl2lph(1:3,1)/norm(Mvxl2lph(:,1));
        dcc = Mvxl2lph(1:3,2)/norm(Mvxl2lph(:,2));
        dcs = Mvxl2lph(1:3,3)/norm(Mvxl2lph(:,3));
        st = norm(Mvxl2lph(:,3));
        PixelSpacing = colvec(sqrt(sum(Mvxl2lph(:,[1 2]).^2)));
        ImageOrientationPatient = [Mvxl2lph(1:3,1)/norm(Mvxl2lph(1:3,1)); Mvxl2lph(1:3,2)/norm(Mvxl2lph(1:3,2))];
        ImagePositionPatient = Mvxl2lph(1:3,:)*[1 1 1 1]';
        for fi = 1:nslices
% Should also force PatientAge and PatentSex
          hdrs{fi}.SliceThickness = st;
          hdrs{fi}.SpacingBetweenSlices = st;
          hdrs{fi}.PixelSpacing = PixelSpacing;
          hdrs{fi}.ImageOrientationPatient = ImageOrientationPatient;
          hdrs{fi}.ImagePositionPatient = ImagePositionPatient + (fi-1)*dcs*st;
          hdrs{fi}.InstanceNumber = fi;
          hdrs{fi}.SliceLocation = (fi-1)*st;
        end

        if isempty(gradwarpinfo)  || ~isfield(gradwarpinfo,'gwtype')
          vol_corr = vol;
        else
          vol_corr = ctx_unwarp_grad(vol,gradwarpinfo.gwtype,gradwarpinfo.unwarpflag,gradwarpinfo.isoctrflag);
        end
        vol_corr_bak = vol_corr; hdrs_bak = hdrs; 

        %  Resample to thinner slices 
        if norm(vol_corr.Mvxl2lph(:,3))>3
          fprintf(1,'Resampling FLAIR acquaision to thinner slices\n');
          [vol_tmp M_tmp] = QD_ctx2mgh(vol_corr); M_tmp = M_RAS_TO_LPH*M_tmp; dims_tmp = size(vol_tmp); dims_tmp2 = dims_tmp; dims_tmp2(3) = dims_tmp(3)*2;
          M_tmp2 = M_tmp; M_tmp2(:,3) = M_tmp(:,3)/2; dr = M_tmp(1:3,:)*[(dims_tmp+1)/2 1]'-M_tmp2(1:3,:)*[(dims_tmp2+1)/2 1]'; M_tmp2(1:3,4) = M_tmp2(1:3,4) - dr;
          vol_ref = ctx_mgh2ctx(zeros(dims_tmp2),M_LPH_TO_RAS*M_tmp2);
          if 0
            vol_flair_t1 = vol_resample(vol_corr,vol_corr_t1,eye(4),2);
            showVol(vol_flair_t1,vol_corr_t1)
            keyboard
          end
          vol_corr = vol_resample(vol_corr_bak,vol_ref,eye(4),2);
          st = norm(vol_corr.Mvxl2lph(:,3)); 
          st = min(st,3); % Fix numeric issue
          if isfield(hdrs{1},'PerFrameFunctionalGroupsSequence') && ~isempty(hdrs{1}.PerFrameFunctionalGroupsSequence.Item_1) % Enhanced dicom?
            if length(hdrs)>1
              fprintf('%s: ERROR: Enhanced dicom should have length(hdrs)==1\n',mfilename);
            end
            fi = 1;
            pos0 = hdrs{fi}.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient;
            for slicenum = 1:nslices
              itemname = sprintf('Item_%d',slicenum);
	      structname = sprintf('hdrs{fi}.PerFrameFunctionalGroupsSequence');
              cmd = sprintf('%s = setfield(%s,''%s'',hdrs{fi}.PerFrameFunctionalGroupsSequence.Item_1);',structname,structname,itemname'); eval(cmd); disp(cmd);
              structname = sprintf('hdrs{fi}.PerFrameFunctionalGroupsSequence.%s.PixelMeasuresSequence.Item_1',itemname);
              cmd = sprintf('%s = setfield(%s,''SliceThickness'',st);',structname,structname'); disp(cmd); eval(cmd); 
              structname = sprintf('hdrs{fi}.PerFrameFunctionalGroupsSequence.%s.PlanePositionSequence.Item_1',itemname);
              pos = pos0; pos(3) = pos0(3)+(slicenum-1)*st;
              cmd = sprintf('%s = setfield(%s,''ImagePositionPatient'',pos);',structname,structname); disp(cmd); eval(cmd); 
            end
          else
            hdrs = repmat({hdrs_bak{1}},[dims_tmp2(3)]);
            for fi = 1:length(hdrs)
              hdrs{fi}.SliceThickness = st; % hdrs{1}.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.SliceThickness in enhanced
              hdrs{fi}.SpacingBetweenSlices = st;
              hdrs{fi}.SliceLocation = (fi-1)*st;
              hdrs{fi}.ImagePositionPatient(3) = hdrs_bak{1}.ImagePositionPatient(3)+(fi-1)*st;
            end
          end
        end

        %        showVol(vol_corr); drawnow;

        % Make sure NQ/LQ app does not perform gradwarp
        for fi = 1:length(hdrs)
          if exist('age') && ~isempty(age), hdr.PatientAge = sprintf('%03dY',age); hdr.PatientSex = sex; end
          if isfield(hdrs{fi},'ImageType')
            hdrs{fi}.ImageType = strrep(hdrs{fi}.ImageType,'ND','DIS3D'); % Should check if Siemens? What about GE or Philips
          end
          if isfield(hdrs{fi},'ManufacturerModelName')
            hdrs{fi}.ManufacturerModelName = 'Dummy';
          end
          if isfield(hdrs{fi},'ManufacturersModelName')
          hdrs{fi}.ManufacturersModelName = 'Dummy';
          end
          if isfield(hdrs{fi},'AccessionNumber') && length(hdrs{fi}.AccessionNumber)>16
            hdrs{fi}.AccessionNumber = hdrs{fi}.AccessionNumber(1:16);
          end
        end

        hdrs_bak = hdrs;

        % Align with re-centered T1 from above
        if isfield(hdrs{1},'PerFrameFunctionalGroupsSequence') && ~isempty(hdrs{1}.PerFrameFunctionalGroupsSequence.Item_1) % Enhanced dicom?
          if ~isfield(hdrs{fi},'PulseSequenceName'), hdrs{fi}.PulseSequenceName = 't2_space_flair'; end
          nslices = size(vol_corr.imgs,3);
          for slicenum = 1:nslices
            itemname = sprintf('Item_%d',slicenum);
            structname = sprintf('hdrs{fi}.PerFrameFunctionalGroupsSequence.%s.PlanePositionSequence.Item_1',itemname);
            pos = eval(sprintf('%s.ImagePositionPatient',structname))+ctr;
            cmd = sprintf('%s = setfield(%s,''ImagePositionPatient'',pos);',structname,structname); disp(cmd); eval(cmd);
          end
        else
          for fi = 1:length(hdrs)
            hdrs{fi}.ImagePositionPatient(:) = hdrs_bak{fi}.ImagePositionPatient(:)-ctr; % This should be hdrs{1}.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1 for enhanced dicom
          end
        end

        if ~exist(outputdir,'dir'), mkdir(outputdir); end

        dirname_flair_out = outputdir;
        dicomwriteVol(vol_corr, dirname_flair_out, hdrs);
        save(fname_flair_info,'hdrs');

        if 0
          vol_t1_dcm = dicomreadVolume_amd(dirname_t1_out);
          vol_flair_dcm = dicomreadVolume_amd(dirname_flair_out); vol_flair_res = vol_resample(vol_flair_dcm,vol_t1_dcm,eye(4),2); showVol(vol_t1_dcm,vol_flair_res); 
        end 

      else
        fprintf(1,'file %s already exists\n',fname_flair_info);
      end
      dirname_flair_out = outputdir;
    end
% Should double-chack that 3D T1 and FLAIR written out are in register
  end
  for phase = 1:2
%  for phase = 2
    cmd = sprintf('singularity run -B %s:/tmp/input1 -B %s:/tmp/output %s NeuroQuant -N -i1 /tmp/input1 --max-measurement-index 1000 --reports NQARA --output-directory /tmp/output',dirname_t1_out,dirname_out,nq4sif);
    fname_done = sprintf('%s/NQ_done.mat',dirname_out); fname_failed = sprintf('%s/NQ_failed.mat',dirname_out);
    if phase==2
      if ~exist('dirname_flair_out','var') || isempty(dirname_flair_out), continue; end
      cmd = sprintf('singularity run -B %s:/tmp/input1 -B %s:/tmp/input2 -B %s:/tmp/output %s LesionQuant -N -i1 /tmp/input1 -i2 /tmp/input2 --max-measurement-index 20 --reports LQBR --output-directory /tmp/output',dirname_t1_out,dirname_flair_out,dirname_out,nq4sif);
      fname_done = sprintf('%s/LQ_done.mat',dirname_out); fname_failed = sprintf('%s/LQ_failed.mat',dirname_out);
    end
    if (exist('sex','var') & exist('age','var')) && (~isempty(sex) && ~isempty(age))
      cmd = sprintf('%s --age %d --sex %s',cmd,age,sex);
    end
    if exist(fname_done,'file')
      fprintf(1,'file %s already exist\n',fname_done);
    else
      [s o e] = jsystem(sprintf('rm %s %s',fname_done,fname_failed));
      disp(cmd)
      [s o e] = jsystem(cmd);
      if s
        disp(s); disp(o); disp(e);
        save(fname_failed,'s','o','e');
      else
        save(fname_done,'s','o','e');
      end
      if ~s
        fprintf(1,'***Done*** %s %s (%0.2f mins)\n',cmd,datestr(now),(now-starttime)*60*24);
      else
        fprintf(1,'***Failed*** %s %s (%0.2f mins)\n',cmd,datestr(now),(now-starttime)*60*24);
      end
    end
  end
  endtime = now;
  %    save(fname_complete,'starttime0','endtime');
catch ME
  fprintf(1,'%s -- Crap!\n',mfilename);
  PrintErrorStruct(ME)
  if isinteractive,keyboard; end
  rethrow(ME)
end

fprintf(1,'%s: ***Completed*** %s (%0.2f mins)\n',mfilename,datestr(now),(now-starttime0)*60*24);

PrintMemoryUsage

% Add code to read in and display results

return

% Display results -- whould read in as ctx structs
[vol_head M] = QD_load_mgh(sprintf('%s/nq/head.mgz',dirname_out));
[vol_seg M] = QD_load_mgh(sprintf('%s/nq/seg.mgz',dirname_out));

% ToDos
%   Make work with nifti files as input
%   Allow for dicom tags to be specified by caller




% Remove done files for ADNI and UKB
flist = dir('/space/syn65/1/data/UKB/NQLQ/*.zip/*/*.mat');
flist = dir('/space/syn65/1/data/ADNI/NQLQ/*/*/*.mat');
for fi = 1:length(flist)
  fname = sprintf('%s/%s',flist(fi).folder,flist(fi).name);
  cmd = sprintf('rm %s',fname);
  fprintf(1,'%s\n',cmd);
  [s o e] = jsystem(cmd);
  if s
    disp([s o e])
  end
end

