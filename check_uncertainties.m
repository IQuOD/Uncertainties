%Look through the netcdf files from NOAA to see if the uncertainty has been
%applied correctly.
clear
pth = '/Volumes/netcdf/';

% %load up the uncertainty csv table
% fid = fopen('iquod_uncertainties_table_WOD_v1.csv');
% c = textscan(fid,'%s%f%s%u32%u32%s%s','delimiter',',','headerlines',1);
% fclose(fid);
% 
% [uinst,ia,ib] = unique(c{3});


%%
%go through each NC file and look at the profile type, year, anything else
%and see if the uncertainties in depth and temperature and salinity have
%been applied.

dirs = dir(pth);
fid = fopen('/Users/cow074/Documents/work_mac/IQuOD/UncertaintyChecks/errors_xbt.txt','a');
for a = 117:length(dirs)
    fnms = dir([pth dirs(a).name '/*xbt*.nc']);
    
    %now go through the netcdf files in each directory
    for b = 1:length(fnms)
        clear t* u* i* z*
        fn = [pth dirs(a).name '/' fnms(b).name]; 
        disp(fn)
        
        %read the country
        country = ncread(fn,'country');
        %uniqueid
        uid = ncread(fn,'wod_unique_cast');
        %date
        dat = ncread(fn,'date');
        %instrument
        t_inst = ncread(fn,'Temperature_Instrument');
        uinst = cellstr(unique(t_inst','rows'));
        
        %read the temperature parts
        t_unc = ncread(fn,'Temperature_uncertainty');
        t_rsize = ncread(fn,'Temperature_row_size');
        if any(t_rsize < 0)
            fprintf(fid,'%s',[fnms(b).name ', Missing values in Temperature_row_size field, file not fully checked.']);
            continue
        end
        %         t_uncal = ncread(fn,'Temperature_uncalibrated');
        
        %         %read the depth parts:
        z = ncread(fn,'z');
        z_unc = ncread(fn,'z_uncertainty');
%         z_rsize = ncread(fn,'z_row_size');
        
%         %         %read the salinity parts
%         try
%             s_unc = ncread(fn,'Salinity_uncertainty');
%             s_rsize = ncread(fn,'Salinity_row_size');
%             s_uncal = ncread(fn,'Salinity_uncalibrated');
%         catch
%         end
        
        %expand the instrument type array:
        inst = cell(size(t_unc));
        uqid = NaN*ones(size(t_unc));
        ind = 1;
        for mm = 1:length(t_rsize)
            indst = ind;
            indend = sum(t_rsize(1:mm));
            rng = indst:indend;
            inst(rng) = cellstr(t_inst(:,mm)');
            uqid(rng) = uid(mm);
            %reset the start point
            ind = indend+1;
        end
        %expand the depth array:

%XBT
        ii = strmatch('XBT',inst);
        if ~isempty(ii)
            %check for Sippican, sub XBT and TSK
            isipp = find(cellfun(@isempty,strfind(inst,'SIPPICAN'))==0);
            usip = find(round(t_unc(isipp)*10) ~=1);
            if ~isempty(usip)
                disp('Sippican Temp XBT uncertainties wrong!')
%                 usip
                uq = unique(uqid(isipp(usip)));
                for mm = 1:length(uq)
                    fprintf(fid,'%s\n',[fnms(b).name ',' num2str(uq(mm)) ',Sippican Temp']);
                end
            end
            itsk = find(cellfun(@isempty,strfind(inst,'TSK'))==0);
            utsk = find(round(t_unc(itsk)*100) ~= 15);
            if ~isempty(utsk)
                disp('TSK XBT Temp uncertainties wrong!')
%                 utsk
                uq = unique(uqid(itsk(utsk)));
                for mm = 1:length(uq)
                    fprintf(fid,'%s\n',[fnms(b).name ',' num2str(uq(mm)) ',TSK Temp']);
                end
            end
            isub = find(cellfun(@isempty,strfind(inst,'SUBMARINE'))==0);
            usub = find(round(t_unc(isub)*100) ~=15);
            if ~isempty(usub)
                disp('Submarine XBT Temp uncertainties wrong!')
%                 usub
                uq = unique(uqid(isub(usub)));
                for mm = 1:length(uq)
                    fprintf(fid,'%s\n',[fnms(b).name ',' num2str(uq(mm)) ',SUBMARINE XBT Temp']);
                end
            end
            %all others should be 0.2
            iall = ones(length(t_unc),1);
            iknown = [isipp; itsk; isub];
            iall(iknown) = 0;
            iothers = find(iall);
            
            uoth = find(round(t_unc(iothers)*10)~=2);
            if ~isempty(uoth)
                disp('Other XBT types Temp uncertainty wrong!')
%                 uoth
                [uq,ia,ib] = unique(uqid(iothers(uoth)));
                for mm = 1:length(uq)
                    fprintf(fid,'%s\n',[fnms(b).name ',' num2str(uq(mm)) ',' ...
                        inst{iothers(ia(mm))}]);
                end
            end
            %             XBT depth
            % we have the idex of all XBT depths (ii):
            if length(z_unc) ~= length(t_unc)
                disp('Depth and temp arrays not same size!')
                keyboard
            end
            i230 = find(z(ii) <= 230);
            id = find(z(ii) > 230);
            u230 = unique(z_unc(ii(i230)));
            if round(u230*10) ~=46
                disp('Depth uncertainty error in <=230m XBTs')
                keyboard
            end
            derr = 0.02*z(ii(id));
            if any(abs(derr - z_unc(ii(id))) > 0.5)
                disp('Depth error in >230m XBTs')
                keyboard
            end
        end
        
        
    end
end
fclose(fid);
