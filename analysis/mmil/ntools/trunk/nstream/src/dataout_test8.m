function dataout_test8

for i = 1:8
    fid = fopen(['dataout.', num2str(i)]);
    offset(i) = fread(fid,1,'int32');
    fclose(fid);
end
offset = min(offset); % left-pad
DataTank = zeros(127,1000);
for Dsp = 1:8
    DspInds = (Dsp-1)*16+[1:16];  % Channel indices
    fid(Dsp) = fopen(['dataout.', num2str(Dsp)]);
    chk = 1;
    while chk
        ts = fread(fid(Dsp),1,'int32');
        if length(ts)
            st = ts-offset+1;
            data = fread(fid(Dsp),[16,928/2/16],'int16');
            DataTank(DspInds,st:st+29-1) = data;
        else
            chk = 0;
        end
    end
end
% 
% >> plot(DataTank(1,1:1e3))
% >> plot(DataTank(1:4,1:1e3))
% >> plot(DataTank(1:4,1:1e3)')
% >> plot(DataTank([1:4 17:21],1:1e3)')
% >> plot(DataTank([1:4 17:20],1:1e3)')
% >> plot(DataTank([1:4 17:20 33:36],1:1e3)')
% >> plot(DataTank([1:4 17:20 33:36],1:1e4)')
% >> plot(DataTank([1:4 17:20 33:36 49:52],1:1e4)')
% >> plot(DataTank([1:4 17:20 33:36 49:52 65:68],1:1e4)')
% >> plot(DataTank([1:4 17:20 33:36 49:52 65:68 81:84],1:1e4)')
% >> plot(DataTank([1:4 17:20 33:36 49:52 65:68 81:84 97:100],1:1e4)')
% >> plot(DataTank([1:4 17:20 33:36 49:52 65:68 81:84 97:100 113:116],1:1e4)')
% >> plot(DataTank([1:4 17:20 33:36 49:52 65:68 81:84 97:100 113:116],1:1e4)')
% 
% 
%     data2 = reshape(permute(data,[2,1,3]),7*16,6e4);
%     for i = 1:size(data2,1)
%         data3(i,:) = data2(i,:)+i*100;
%     end
%     plot(data3')
