%Look through the netcdf files from NOAA to see if the uncertainty has been
%applied correctly.
clear
pth = '/Volumes/netcdf/';
% pth = '/Users/cow074/Documents/work_mac/IQuOD/UncertaintyChecks/data/';
%load up the uncertainty csv table
fid = fopen('iquod_uncertainties_table_WOD_v1.csv');
c = textscan(fid,'%s%f%s%u32%u32%s%s','delimiter',',','headerlines',1);
fclose(fid);


[allinst,~,iali] = unique(c{3});
var = c{1};
unc = c{2};
start = c{4};
modf = c{6};

%%
%go through each NC file and look at the profile type, year, anything else
%and see if the uncertainties in depth and temperature and salinity have
%been applied.

dirs = dir(pth);
fid = fopen('/Users/cow074/Dropbox/UncertaintyChecks/errors_osd_ctd_19Oct.txt','a');
for a = 103:length(dirs)
    fnm = dir([pth dirs(a).name '/*osd*.nc']);
    fnm2 = dir([pth dirs(a).name '/*ctd*.nc']);
    fnms = [fnm;fnm2];
    
    %now go through the netcdf files in each directory
    for b = 1:length(fnms)
        clear t_inst uinsf ibk calib* psal* issal
        fn = [pth dirs(a).name '/' fnms(b).name];
        disp(fnms(b).name)
        
        %read the country
        try
            country = ncread(fn,'country');
        catch
            fprintf(fid,'%s\n',['File corrupt?: ' fn]);
            continue
        end
        %uniqueid
        uid = ncread(fn,'wod_unique_cast');
        %date
        dat = ncread(fn,'date');
        %instrument
        try
            t_inst = ncread(fn,'Temperature_Instrument');
            [uinstf,~,ibk] = unique(t_inst','rows');
            uinstf = cellstr(uinstf);
        catch
            fprintf(fid,'%s\n',[fnms(b).name ' has no Temperature_Instrument variable']);
            uinstf = cellstr('');
            ibk = ones(size(dat));
        end
        
        %read the temperature parts
        t = ncread(fn,'Temperature');
        t_unc = ncread(fn,'Temperature_uncertainty');
        t_rsize = ncread(fn,'Temperature_row_size');
        try
            calibt = ncread(fn,'Temperature_uncalibrated');
        catch
            %1=uncalibrated
            calibt = ones(size(t));
        end
        
        if isnan(sum(t_rsize))
            disp('NaN in t_row_size')
            t_rsize(isnan(t_rsize)) = 0;
        end
        if sum(t_rsize) ~= length(t_unc)
            disp('Temp row size and t_unc do not match')
            continue
        end
        
        %         %read the depth parts:
        z = ncread(fn,'z');
        z_unc = ncread(fn,'z_uncertainty');
        z_rsize = ncread(fn,'z_row_size');
        try
            calibz = ncread(fn,'Pressure_uncalibrated');
        catch
            %1=uncalibrated
            calibz = ones(size(z));
        end
        %         %         %read the salinity parts
        try
            psal = ncread(fn,'Salinity');
            issal = 1;
        catch
            psal = [];
            psal_rsize = [];
            psal_unc = [];
            issal = 0;
        end
        try
            psal_rsize = ncread(fn,'Salinity_row_size');
            psal_unc = ncread(fn,'Salinity_uncertainty');
        catch
            disp('missing salinity info in file')
            psal_unc = NaN*ones(size(psal));
        end
        try
            calibps = ncread(fn,'Salinity_uncalibrated');
        catch
            %1=uncalibrated
            calibps = ones(size(psal));
        end
        
        if sum(t_rsize) ~= length(t_unc)
            disp('Temp row size and t_unc do not match')
            continue
        end
        
        if isnan(sum(z_rsize))
            disp('NaN in z_row_size')
            z_rsize(isnan(z_rsize)) = 0;
        end
        if sum(z_rsize) ~= length(z_unc)
            disp('Depth row size and z_unc do not match')
            continue
        end
        
        %figure out start/end indices of each profile:
        indst = [1;cumsum(double(t_rsize(1:end-1)))+1];
        indet = cumsum(double(t_rsize));
        indsps = [1;cumsum(double(psal_rsize(1:end-1)))+1];
        indeps = cumsum(double(psal_rsize));
        indsz = [1;cumsum(double(z_rsize(1:end-1)))+1];
        indez = cumsum(double(z_rsize));
        
        %let's see what we can match:
        iknown = 9*ones(size(uinstf));%OSD default
        for mm = 1:length(allinst)
            iin = find(cellfun(@isempty,strfind(uinstf,upper(allinst{mm})))==0);
            %flag this as known.
            iknown(iin) = mm;
        end
        
        %cycle through each data type
        for mm = 1:length(iknown)
            %             if iknown(mm) == 0
            %                 %no values assigned, what is this datatype?
            %                 disp('Unknown data type in this file')
            %                 uq = uid(ibk == mm);
            %                 for nn = 1:length(uq)
            %                     fprintf(fid,'%s\n',[fnms(b).name ',' num2str(uq(nn)) ',' ...
            %                         uinstf{mm} ',No data type']);
            %                 end
            %             else %check the known types
            
            %check for temperature
            iunct = find(iali == iknown(mm));
            iin = find(cellfun(@isempty,strfind(var(iunct),'temperature')) == 0);
            iis = find(cellfun(@isempty,strfind(var(iunct),'salinity')) == 0);
            iiz = find(cellfun(@isempty,strfind(var(iunct),'depth')) == 0);
            uq = uid(ibk == mm);
            %MBT
            if ~isempty(strfind(uinstf{mm},'MBT'))
                keyboard
            end
            
            %CTD
            if ~isempty(strmatch('CTD',uinstf{mm})) |...
                    ~isempty(strmatch('XCTD',uinstf{mm}))
                rng = [];uuid=[];datn=[];
                ist = indst(ibk==mm); ine = indet(ibk==mm);
                %delete zero indices.
                ist(t_rsize(ibk==mm) == 0) = NaN;
                ine(t_rsize(ibk==mm) == 0) = NaN;
                d = dat(ibk==mm);
                if isempty(strfind(uinstf{mm},'XCTD'))
                    calib = [];
                    cal = calibt(ibk == mm);
                end
                for nn = 1:length(ist)
                    if isnan(ist(nn))
                        continue
                    end
                    rng = [rng;(ist(nn):ine(nn))'];
                    uuid = [uuid;repmat(uq(nn),ine(nn)-ist(nn)+1,1)];
                    datn = [datn;repmat(d(nn),ine(nn)-ist(nn)+1,1)];
                    if isempty(strfind(uinstf{mm},'XCTD'))
                        calib = [calib;repmat(cal(nn),ine(nn)-ist(nn)+1,1)];
                    end
                end
                %figure out which category of temperature error:
                %1 = uncalibrated

                if  ~isempty(strfind(uinstf{mm},'XCTD'))
                    cttype = {'ictd','i1998'};
                    i1998 = datenum(num2str(datn),'yyyymmdd') > datenum('19971231','yyyymmdd');
                    ictd = ones(size(datn));
                    ictd(i1998) = 0;
                    i1998 = find(i1998);ictd = find(ictd);
                else
                    icalib = calib ~= 1;
                    i1980 = datenum(num2str(datn),'yyyymmdd') > datenum('19791231','yyyymmdd');
                    ictd = ones(size(calib));
                    ictd(icalib) = 0; ictd(i1980) = 0;
                    icalib = find(icalib);i1980 = find(i1980);ictd = find(ictd);
                    cttype = {'ictd','i1980','icalib'};
                end
                
                
                %temperature
                for ictd_type = 1:length(cttype)
                    eval(['ii = ' cttype{ictd_type} ';'])
                    if ~isempty(ii)
                        utm = find((round(t_unc(rng(ii))*100) ~= round(unc(iunct(iin(ictd_type)))*100)) ...
                            & ~isnan(t(rng(ii))));
                        if ~isempty(utm)
                            disp([uinstf{mm} ' Temp uncertainties wrong!'])
                            [printonce,ip] = unique(uuid(utm),'stable');
                            for nn = 1:length(printonce)
                                fprintf(fid,'%s\n',[fnms(b).name ',' ...
                                    num2str(printonce(nn)) ',' uinstf{mm} ',' ...
                                    num2str(t_unc(rng(ii(utm(ip(nn)))))) ',Temp_unc']);
                            end
                        end
                    end
                end
                %salinity
                if issal
                    rng = [];uuid=[];
                    ist = indsps(ibk==mm); ine = indeps(ibk==mm);
                    %delete zero indices.
                    ist(psal_rsize(ibk==mm) == 0) = NaN;
                    ine(psal_rsize(ibk==mm) == 0) = NaN;
                    if isempty(strfind(uinstf{mm},'XCTD'))
                        cal = calibps(ibk == mm);
                        calib = [];
                    end
                    for nn = 1:length(ist)
                    if isnan(ist(nn))
                        continue
                    end
                        rng = [rng;(ist(nn):ine(nn))'];
                        uuid = [uuid;repmat(uq(nn),ine(nn)-ist(nn)+1,1)];
                        if isempty(strfind(uinstf{mm},'XCTD'))
                            calib = [calib;repmat(cal(nn),ine(nn)-ist(nn)+1,1)];
                        end
                    end
                    
                    if ~isempty(strfind(uinstf{mm},'XCTD'))
                        i1998 = datenum(num2str(datn),'yyyymmdd') > datenum('19971231','yyyymmdd');
                        ictd = ones(size(datn));
                        ictd(i1998) = 0;
                        ictd = find(ictd);i1998 = find(i1998);
                        cttype = {'ictd','i1998'};
                     else
                       icalib = calib ~= 1;
                        ictd = ones(size(calib));
                        ictd(icalib) = 0;
                        icalib = find(icalib);ictd = find(ictd);
                        cttype = {'ictd','icalib'};
                    end
                    for ictd_type = 1:2
                        eval(['ii = ' cttype{ictd_type} ';'])
                        if isempty(rng)
                            continue
                        end
                        if ~isempty(ii)
                            utm = find((round(psal_unc(rng(ii))*100) ~= round(unc(iunct(iis(ictd_type)))*100)) ...
                                & ~isnan(psal(rng(ii))));
                            if ~isempty(utm)
                                %                                     keyboard
                                disp([uinstf{mm} ' PSAL uncertainties wrong!'])
                                [printonce,ip] = unique(uuid(utm),'stable');
                                for nn = 1:length(printonce)
                                    fprintf(fid,'%s\n',[fnms(b).name ',' ...
                                        num2str(printonce(nn)) ',' uinstf{mm} ',' ...
                                        num2str(psal_unc(rng(ii(utm(ip(nn)))))) ',Psal_unc']);
                                end
                            end
                        end
                    end
                end
                
                %depth
                rng = [];uuid=[];datn=[];
                ist = indsz(ibk==mm); ine = indez(ibk==mm);d = dat(ibk==mm);
                %delete zero indices.
                ist(z_rsize(ibk==mm) == 0) = NaN;
                ine(z_rsize(ibk==mm) == 0) = NaN;
                if isempty(strfind(uinstf{mm},'XCTD'))
                    calib = [];
                    cal = calibz(ibk==mm);
                end
                for nn = 1:length(ist)
                    if isnan(ist(nn))
                        continue
                    end
                    rng = [rng;(ist(nn):ine(nn))'];
                    uuid = [uuid;repmat(uq(nn),ine(nn)-ist(nn)+1,1)];
                    datn = [datn;repmat(d(nn),ine(nn)-ist(nn)+1,1)];
                    if isempty(strfind(uinstf{mm},'XCTD'))
                        calib = [calib;repmat(cal(nn),ine(nn)-ist(nn)+1,1)];
                    end
                end
                
                if isempty(strfind(uinstf{mm},'XCTD'))
                    icalib = calib ~= 1;
                    ictd = ones(size(calib));
                    ictd(icalib) = 0;
                    i1980 = datenum(num2str(datn),'yyyymmdd') > datenum('19791231','yyyymmdd');
                    ictd(i1980) = 0;
                    icalib = find(icalib);ictd = find(ictd);i1980 = find(i1980);
                    icalib = unique([i1980;icalib]); %both have the same uncertainty in depth
                    cttype = {'ictd','icalib'};
                else
                    ictd = ones(size(datn));
                    i1998 = datenum(num2str(datn),'yyyymmdd') > datenum('19971231','yyyymmdd');
                    ictd(i1998) = 0;
                    ictd = find(ictd);i1998 = find(i1998);
                    cttype = {'ictd','i1998'};
                end
                for ictd_type = 1:2
                        if isempty(rng)
                            continue
                        end
                    eval(['ii = ' cttype{ictd_type} ';'])
                    if ~isempty(ii)
                        utm = find(abs(round(z_unc(rng(ii))*1000) - round(z(rng(ii))*unc(iunct(iiz(ictd_type)))*10)) >1 ...
                            & ~isnan(z(rng(ii))));
                        if ~isempty(utm)
                            %                                     keyboard
                            disp([uinstf{mm} ' Depth uncertainties wrong!'])
                            [printonce,ip] = unique(uuid(utm),'stable');
                            for nn = 1:length(printonce)
                                fprintf(fid,'%s\n',[fnms(b).name ',' ...
                                    num2str(printonce(nn)) ',' uinstf{mm} ',' ...
                                    num2str(z_unc(rng(ii(utm(ip(nn)))))) ',z_unc']);
                            end
                        end
                    end
                end
                
                
            else
                %find the uncertainty and compare
                rng = [];uuid=[];
                ist = indst(ibk==mm); ine = indet(ibk==mm);
                %delete zero indices.
                ist(t_rsize(ibk==mm) == 0) = NaN;
                ine(t_rsize(ibk==mm) == 0) = NaN;
                for nn = 1:length(ist)
                    if isnan(ist(nn))
                        continue
                    end
                    rng = [rng;(ist(nn):ine(nn))'];
                    uuid = [uuid;repmat(uq(nn),ine(nn)-ist(nn)+1,1)];
                end
                %temperature
                utm = find((round(t_unc(rng)*100) ~= round(unc(iunct(iin))*100)) ...
                    & ~isnan(t(rng)));
                if ~isempty(utm)
                    disp([uinstf{mm} ' Temp uncertainties wrong!'])
                    [printonce,ip] = unique(uuid(utm),'stable');
                    for nn = 1:length(printonce)
                        fprintf(fid,'%s\n',[fnms(b).name ',' ...
                            num2str(printonce(nn)) ',' uinstf{mm} ',' ...
                            num2str(t_unc(rng(utm(ip(nn))))) ',Temp_unc']);
                    end
                end
                %depth - percentage or m?
                %assume Tim has interpreted it as m
                %                     derr = unc(iunct(iiz))*z(rng)/100;
                rng = [];uuid=[];
                ist = indsz(ibk==mm); ine = indez(ibk==mm);
                %delete zero indices.
                ist(z_rsize(ibk==mm) == 0) = NaN;
                ine(z_rsize(ibk==mm) == 0) = NaN;
                for nn = 1:length(ist)
                    if isnan(ist(nn))
                        continue
                    end
                    rng = [rng;(ist(nn):ine(nn))'];
                    uuid = [uuid;repmat(uq(nn),ine(nn)-ist(nn)+1,1)];
                end
                derr = unc(iunct(iiz)) * z(rng)/100;
                utz = find(round(derr*1000) ~= round(z_unc(rng)*1000)); % & abs(derr - z_unc(rng)) < 99999);
                if ~isempty(utz)
                    disp(['Depth error in ' uinstf{mm}])
                    [printonce,ip] = unique(uuid(utz),'stable');
                    for nn = 1:length(printonce)
                        fprintf(fid,'%s\n',[fnms(b).name ',' ...
                            num2str(printonce(nn)) ',' uinstf{mm} ',' ...
                            num2str(z_unc(rng(utz(ip(nn))))) ',z_unc']);
                    end
                end
            end
            
        end
    end
end
    fclose(fid);
