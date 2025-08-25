function dicomwriteVol(vol, outputdir, hdrs, varargin) 

% vol - volume (as ctx struct) to be written into dicom
% outputdir - directory to save created dicoms
% hdrs -- metadata structs (from dicominfo) -- should have same length and ordering as vol slices
% should also allow for name/value pairs in varargin

SeriesInstanceUID = dicomuid;

dat = vol.imgs;

fprintf('%s -- %s | Writing DICOM images...\n', datestr(now), mfilename);

if length(hdrs)==1 && vol.dimd>1 % Enhanced disom?
  i = 1;
  metadata = hdrs{i};
  metadata.SeriesInstanceUID = SeriesInstanceUID;
  metadata.SOPInstanceUID = dicomuid;
  fname_out = sprintf('%s/im%3.4i.dcm',outputdir,i);
  switch metadata.BitDepth
    case 8
      metadata.BitsAllocated = 8;
      metadata.BitsStored = 8;
      metadata.HighBit = 7;
      dicomwrite(uint8(permute(squeeze(dat),[2 1 4 3])), fname_out, metadata, 'CreateMode', 'copy' ); % Not sure if writing private tags works with NQ
    case 16
      metadata.BitsAllocated = 16;
      metadata.BitsStored = 16;
      metadata.HighBit = 15;
      dicomwrite(uint16(permute(squeeze(dat),[2 1 4 3])), fname_out, metadata, 'CreateMode', 'copy' );
    otherwise % Need some way to handle metadata.BitsStored == 12 , check into difference between metadata.BitsAllocated and metadata.BitsStored
      metadata.BitsAllocated = 16;
      metadata.BitsStored = 16;
      metadata.HighBit = 15;
      dicomwrite(uint16(permute(squeeze(dat),[2 1 4 3])), fname_out, metadata, 'CreateMode', 'copy' );
  end
%  vol_in = permute(squeeze(dicomread(fname_out)),[2 1 3 4]);;
  fprintf('file %s written (%s)\n',fname_out,datestr(now));
else
  for i = 1:vol.dimd
    metadata = hdrs{i};
    metadata.SeriesInstanceUID = SeriesInstanceUID;
    metadata.SOPInstanceUID = dicomuid;
    fname_out = sprintf('%s/im%3.4i.dcm',outputdir,i);
    switch metadata.BitDepth
      case 8
        metadata.BitsAllocated = 8;
        metadata.BitsStored = 8;
        metadata.HighBit = 7;
        dicomwrite( uint8( permute(squeeze(dat(:,:,i,:)), [2 1 3]) ), fname_out, metadata, 'CreateMode', 'copy' ); % Not sure if writing private tags works with NQ
      case 16
        metadata.BitsAllocated = 16;
        metadata.BitsStored = 16;
        metadata.HighBit = 15;
        dicomwrite( uint16( permute(squeeze(dat(:,:,i,:)), [2 1 3]) ), fname_out, metadata, 'CreateMode', 'copy' );
      otherwise
        metadata.BitsAllocated = 16;
        metadata.BitsStored = 16;
        metadata.HighBit = 15;
        dicomwrite( uint16( permute(squeeze(dat(:,:,i,:)), [2 1 3]) ), fname_out, metadata, 'CreateMode', 'copy' );
    end
    fprintf('file %s written (%s)\n',fname_out,datestr(now));
  end
end
fprintf('\n%s - DICOM directory %s written (%d files)\n',mfilename,outputdir,vol.dimd);

