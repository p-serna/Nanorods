function Movie=ReadLongMovie(filename, Nframes)

% Use low-level File I/O to read the file
info = imfinfo(filename);
info = info(1, :);
fp = fopen(filename , 'rb');
% The StripOffsets field provides the offset to the first strip. Based on
% the INFO for this file, each image consists of 1 strip.
fseek(fp, info.StripOffsets, 'bof');
% Assume that the image is 16-bit per pixel and is stored in big-endian format.
% Also assume that the images are stored one after the other.

% Read the movie
Movie=zeros(info.Height,info.Width,Nframes);
%Movie=zeros(info.Height,info.Width,Nframes,'uint16');
for idx=1:Nframes
    Movie(:,:,idx)=fread(fp, [info.Width info.Height], 'uint16', 0, 'ieee-be')';
end
fclose(fp);

