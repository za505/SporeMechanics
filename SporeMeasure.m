%%SporeMeasure.m
%Author: Zarina Akbary
%Date: 30 June 2022
%Purpose: to measure the long and short axis of B. subtilis PY79 spores in
%isosmotic and hyperosmotic media

clear, close all

%User Input
dirname = ['/Users/zarina/Downloads/EJProject'];
basename = {'20220627_WT_turgid', '20220627_WT_hypershock', '20220707_cotXYZ_turgid', '20220707_cotXYZ_hypershock', '20220707_cotE_turgid', '20220707_cotE_hypershock'};
labels = {'WT isosmotic', 'WT hypershock', 'cotE isosmotic', 'cotE hypershock', 'cotXYZ isosmotic', 'cotXYZ hypershock'};

%okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [204, 121, 167], [213, 94, 0], [0, 114, 178], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);
okabeIto = [okabeIto, okabeIto];

sm = 3; %parameter for edge detection
dr=1;%Radius of dilation before watershed (water leaves cell) 
lscale = 0.08; %Microns per pixel
minA = 15; %minimum area for initial filter
maxA = 70; %maximum area for final filter
minF = 50; %minimum area for final filter
maxF = 250;%maximum area for final filter
troubleshooting = 0; %visualize cells

%% load images and measure spores
[boundaries, centroid, lCell, nCells, pixloc, wCell] = sporeMeasure(dirname, basename, dr, lscale, maxA, maxF, minA, minF, sm, troubleshooting);

%% plot the dimensions in scatter plot
p = cell(length(basename),1);

figure, hold on
for i = 1:length(basename) %bc using an array
    [nrow, ncol] = size(lCell{i});
    x = lCell{i};
    y = wCell{i};
    p{i} = scatter(lCell{i}, wCell{i}, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i})
end
xlabel('axis 1 (\mum)')
ylabel('axis 2 (\mum)')
legend([p{1}(1), p{2}(1), p{3}(1), p{4}(1), p{5}(1), p{6}(1)], labels, 'Location', 'southeast')

%% plot the dimensions as subplots
p = cell(length(basename),1);
j = 1;
strain = {'', 'WT', '', 'cotE', '', 'cotXYZ'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 20), hold on
for i = 1:length(basename) %bc using an array
    subplot(1, 3, j), hold on
    [nrow, ncol] = size(lCell{i});
    x = lCell{i};
    y = wCell{i};
    p{i} = scatter(lCell{i}, wCell{i}, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i})
    xlabel('axis 1 (\mum)')
    ylabel('axis 2 (\mum)')
  
    if mod(i,2)==0
        j = j+1;
        legend([p{i-1}(1), p{i}(1)], {'Isosmotic', 'Hyperosmotic'}, 'Location', 'southeast')
        title(strain{i})
    end
end

%% Determine covariance between axis 1 and axis 2
p = cell(length(basename),1);

figure, hold on
subplot(1,2,1), hold on
for i = 1:length(basename) %bc using an array
    [nrow, ncol] = size(lCell{i});
    x1 = reshape(lCell{i}, [nrow*ncol, 1]);
    p{i} = histogram(x1, 'EdgeColor', okabeIto{i}, 'FaceColor', okabeIto{i}, 'BinWidth', 0.2)
end
xlim([0.2, Inf])
title('Axis 1 (\mum)')
legend([p{1}(1), p{2}(1), p{3}(1), p{4}(1), p{5}(1), p{6}(1)], labels, 'Location', 'southeast')

subplot(1,2,2), hold on
for i = 1:length(basename) %bc using an array
    [nrow, ncol] = size(lCell{i});
    x2 = reshape(wCell{i}, [nrow*ncol, 1]);
    p{i} = histogram(x2, 'EdgeColor', okabeIto{i}, 'FaceColor', okabeIto{i}, 'BinWidth', 0.2)
end
xlim([0.2, Inf])
title('Axis 2 (\mum)')
legend([p{1}(1), p{2}(1), p{3}(1), p{4}(1), p{5}(1), p{6}(1)], labels, 'Location', 'southeast')

%% Functions
function [boundaries, centroid, lCell, nCells, pixloc, wCell] = sporeMeasure(dirname, basename, dr, lscale, maxA, maxF, minA, minF, sm, troubleshooting)

    srow = length(basename);
    boundaries = cell(srow, 1);
    centroid = cell(srow, 1);
    lCell = cell(srow, 1);
    nCells = cell(srow,1);
    pixloc = cell(srow,1);
    wCell = cell(srow, 1);

    for b = 1:srow

        base = basename{b};
        cd([dirname '/' base]);
        directory = dir('*.tif'); %making directory of tif files
        N = length(directory);

        ncells = nan(N, 1); 
        centroids = cell(N,1);
        pixels = cell(N,1);
        B = cell(N,1);
        pcell = nan(N, 1);
        DS = nan(N, 1);
        mlines = cell(N,1);
        acell = nan(N, 1);
        lcell = nan(N, 1);
        wcell = nan(N, 1);

    for n=1:N

        imname = directory(n).name; %changing tiff file name 
        im = imread(imname);
        [imM, imN] = size(im);

        %figure(1), imshow(im, [])
        %figure(2), imhist(im), pause, close

        %De-speckle image
        im=medfilt2(im);

        %Normalize images
        ppix=0.5; %default 0.5
        im=norm16bit(im,ppix);

        %figure(1), imshow(im)
        %figure(2), imhist(im), pause, close

        %enhance contrast
        imc = im; %imc = de-speckled image

        [imcounts,bins]=imhist(imc); 
        [imcounts,idx]=sort(imcounts);
        bins=bins(idx);
        thresh1=bins(end-1);
        if thresh1==65535
            thresh1=bins(end);
        end

        imc=imadjust(imc,[thresh1/65535 1],[]);   

        %Find edges
        [ed2,thresh2]=edge(imc,'canny',[],sm*sqrt(2));

%         [y, x] = find(ed2==1);
%         figure(1)
%         imshow(imc), hold on, plot(x,y,'.r'), pause, close

        %Clean image
        cc=bwconncomp(ed2,8); %default 8
        stats=regionprops(cc,imc,'Area','MeanIntensity');
        idx=find([stats.Area]>minA&[stats.Area]<maxA);
        ed2=ismember(labelmatrix(cc),idx);

        %Close gaps in edges
        despurred=bwmorph(ed2,'spur');
        spurs=ed2-despurred;
        [spy,spx]=find(spurs);
        for k=1:length(spx)
            ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1)=ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1)+rot90(ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1),2);
            ed2(spy(k),spx(k))=1;
        end
        ed2=bwmorph(ed2,'bridge'); 

        se=strel('disk',dr);
        ed2=imdilate(ed2,se);
        ed2=bwmorph(ed2,'thin',2);

        %Identify cells based on size and intensity
        ed3=~ed2;
        ed3(1,:)=ones(1,imN);
        ed3(end,:)=ones(1,imN);
        ed3(:,1)=ones(imM,1);
        ed3(:,end)=ones(imM,1);

        cc=bwconncomp(ed3,4);
        stats=regionprops(cc,imc,'Area','MeanIntensity');

        idx=find([stats.Area]>minF&[stats.Area]<maxF);
        ed4=ismember(labelmatrix(cc),idx);

        %Find cell areas and centroids
        bw=bwmorph(ed4,'thicken');
        [extrace,bw]=bwboundaries(bw,4,'noholes');
        stats=regionprops(bw, 'Area', 'Centroid','PixelIdxList');

        %figure(1), imshow(bw), pause, close

        ncells(n,1)=length(extrace);

        %save pixel id list
        for k=1:ncells(n)
            acell(n,k)=stats(k).Area;
            centroids{n,k}=stats(k).Centroid;
            pixels{n,k}=stats(k).PixelIdxList;

            P = extrace(k);

            rP=[P{1}(:,2),P{1}(:,1)];
            px=[rP(1:end-1,1);rP(1:end-1,1);rP(:,1)];
            py=[rP(1:end-1,2);rP(1:end-1,2);rP(:,2)];
            sp=length(rP);
            dS=sqrt(diff(px).^2+diff(py).^2);
            S=[0 cumsum(dS)'];

            px=csaps(S,px,0.05,S);
            py=csaps(S,py,0.05,S);

            px=px(sp+1:2*sp);
            py=py(sp+1:2*sp);

            px(end)=px(1);
            py(end)=py(1);

            dS=sqrt(diff(px).^2+diff(py).^2);
            S=[0 cumsum(dS)];
            ls=length(S);
            DS(n,k)=S(end)/(ls-1);
            Sn=(0:DS(n,k):S(end));
            nx=spline(S,px,Sn);
            ny=spline(S,py,Sn);

            B{n,k}=[nx',ny'];

            X=B{n,k}(:,1);
            Y=B{n,k}(:,2);   

            [sX,~]=size(X);

            %Find poles
            [X,Y,pcell(n,k)]=polefinder(X,Y);

            %Create mesh
            npts=min(pcell(n,k),sX-pcell(n,k)+1);
            S=(0:DS(n,k):(sX-1)*DS(n,k));

            s1=(0:S(pcell(n,k))/(npts-1):S(pcell(n,k)));
            s2=(S(pcell(n,k)):(S(end)-S(pcell(n,k)))/(npts-1):S(end));
            xc1=spline(S(1:pcell(n,k)),X(1:pcell(n,k)),s1);
            xc2=spline(S(pcell(n,k):end),X(pcell(n,k):end),s2);
            yc1=spline(S(1:pcell(n,k)),Y(1:pcell(n,k)),s1);
            yc2=spline(S(pcell(n,k):end),Y(pcell(n,k):end),s2);
            xc2=fliplr(xc2);
            yc2=fliplr(yc2);

            %Calculate midline
            mlines{n,k}=[(xc1+xc2)'/2,(yc1+yc2)'/2];
            dxy=diff(mlines{n,k}).^2;
            dl=sqrt(dxy(:,1)+dxy(:,2));
            lcell(n,k)=sum(dl);

            %Calculate width
            ls=[0 cumsum(dl)'];
            [~,mpos1]=min(abs(ls/lcell(n,k)-0.25));
            [~,mpos2]=min(abs(ls/lcell(n,k)-0.75));

            widths=sqrt((xc1-xc2).^2+(yc1-yc2).^2);
            wcell(n,k)=(widths(mpos1)+widths(mpos2))/2;

        end

        %Dimsionalize the variables
        lcell(n,:)=lcell(n, :)*lscale;
        wcell(n,:)=wcell(n, :)*lscale;
        acell(n,:)=acell(n, :)*lscale^2;

        lcell(lcell==0)=NaN;
        wcell(wcell==0)=NaN;
        acell(acell==0)=NaN;
        pcell(pcell==0)=NaN;

        %visualize the tracking
        if troubleshooting==1
            figure, imshow(im), hold on       
            for k=1:ncells(n)
                plot(B{n,k}(:,1),B{n,k}(:,2),'-r')
            end
            pause, close all
        elseif troubleshooting==2
            for k=1:ncells(n)
                figure, imshow(im), hold on
                plot(B{n,k}(:,1),B{n,k}(:,2),'-r')
                pause, close all
            end           
        end

    end

    boundaries{b} = B;
    centroid{b} = centroids;
    lCell{b} = lcell;
    nCells{b} = ncells;
    pixloc{b} = pixels;
    wCell{b} = wcell;

    end
end