%  ntools_procCleanDecIEEG processes decieeg file to give line-noise-filtered files
%
%  Inputs:  	recording_filename_root  Recording prefix or numbers
%
%
function ntools_procCleanDecIEEG(recording_filename_root)

    global experiment

    recording_filename_root = make_filename_root(recording_filename_root); %removes suffix if necess

    FS = experiment.processing.ieeg.sample_rate;
    CH = length(experiment.channels);
    clear fks
    for i = 1:length(experiment.processing.ieeg.linefilter.frequencies)
        fks(i,:) = experiment.processing.ieeg.linefilter.frequencies(i) + [-10,10];
    end

    tic
    T = 0.7;

    if isfile([recording_filename_root '.decieeg.dat'])

        ieeg_fid = fopen([recording_filename_root '.decieeg.dat'], 'r');
        clean_fid = fopen([recording_filename_root '.cleandecieeg.dat'], 'w');

        chk=1; N=0;
        while(chk)
            tic
            N = N+1;
            disp(['Clean IEEG: Loop ' num2str(N)]);
            data = fread(ieeg_fid,[CH,FS*T],'short=>single');
            clean = zeros(CH,size(data,2),'single');
            if(size(data,2))
                for ii = 1:size(data,1)
                    cl = linefilter(data(ii,:),[min(T,size(data,2)./FS),10],FS,fks);
                    clean(ii,1:size(data,2)) = cl;
                end
                fwrite(clean_fid,clean,'float');
            end
            if (size(data,2) < FS*T); chk = 0; end
            toc
        end
        fclose(ieeg_fid);
        fclose(clean_fid);
    else
        disp('No decieeg file found.  Skipping.')
    end
end
