clear all
%---------------------------------------------------------------------------------------

%MATLAB code used to find out nulcear/perinuclear distribution in the paper
%"Matrix stiffening promotes perinuclear clustering of mitochondria"
%published in MBoC 2024. Written by Basil Thurakkal.
%---------------------------------------------------------------------------------------


%This code rquires the user to input 3 binary tiff files namely, Cell oultine,
%Nucleus and mitochodrial images. The respective organelles should be white while background value is zero. 
%mitochondria binary input file - background should be 255 and signal should be 0.
%outline binary input file -background should be 0 and signal should be 255. 
%nucleus binary input file - background should be 0 and signal should be 255

Outline_Dir = uigetdir; %Cell outline directory
Outline_Files = dir(fullfile(Outline_Dir,'*.tif'));
Nucleus_Dir = uigetdir; %Nucleus directory
Nucleus_Files = dir(fullfile(Nucleus_Dir,'*.tif'));
Mito_Dir = uigetdir; %Mito directory
Mito_Files = dir(fullfile(Mito_Dir,'*.tif'));

for x = 1:length(Outline_Files)

    FileName = Outline_Files(x).name;
    fullFileName = fullfile(Outline_Dir, FileName);
    fprintf(1, 'Now reading %s\n', FileName);
    cell_outline = imread(fullFileName);
    
    Nucl_file_name = Nucleus_Files(x).name;
    fullname_nucl = fullfile(Nucleus_Dir,Nucl_file_name);
    nucleus = imread(fullname_nucl);

    Mito_file_name = Mito_Files(x).name;
    fullname_mito = fullfile(Mito_Dir,Mito_file_name);
    mito = imread(fullname_mito);
    
    %%
    cell_outline(cell_outline>1) = 1; 
    nucleus(nucleus>1) = 1;
    mito(mito>1) = 1;
    %cell_outline = 1-cell_outline;
    mito = 1-mito;
    mito = mito.*cell_outline;
    %%
    nucleus = 1-nucleus;
    %%
    imagesc(cell_outline);
    %%
    imagesc(nucleus);
    %%
    imagesc(mito);
    %%
    cell_outline_nucl= cell_outline.*nucleus;
    imagesc(cell_outline_nucl);
    %%
    nucleus = 1- nucleus;
    CC_nucleus = bwconncomp(nucleus);
    Centroid_nucleus = regionprops(CC_nucleus,'centroid');
    %%
    centroid_floor = floor(cell2mat(struct2cell(Centroid_nucleus)));
    %% Cell outlne distribution. Cell outline/shape distribution is calculated to use as a benchmark for mitochondrial distribution
     outline_matrix = [];
     for i = 1:size(nucleus,1)
         for j = 1:size(nucleus,2)
             if cell_outline_nucl(i,j)==1
                outline_tmp = sqrt(((i-centroid_floor(2))^2)+((j-centroid_floor(1))^2)); %findng distance for each pixels from the nucleus centroid
                outline_matrix = [outline_matrix outline_tmp];
                i
             else
                 continue
             end
         end
     end
     %% Mitochondrial distribution
      mito_matrix = [];
     for i = 1:size(nucleus,1)
         for j = 1:size(nucleus,2)
             if mito(i,j)==1
                mito_tmp = sqrt(((i-centroid_floor(2))^2)+((j-centroid_floor(1))^2)); %findng distance for each mitochondrialpixels from the nucleus centroid
                mito_matrix = [mito_matrix mito_tmp];
                i
             else
                 continue
             end
         end
     end
     %% Normalization
    [N,edges] = histcounts(mito_matrix, 'Normalization','pdf');
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    
    
    [N1,edges1] = histcounts(outline_matrix, 'Normalization','pdf');
    edges1 = edges1(2:end) - (edges1(2)-edges1(1))/2;
    %%
    figure()
    plot(edges, N);
    hold on
    plot(edges1, N1);
    %%

    csvwrite(strrep(fullFileName,'.tif','_mito.csv'),mito_matrix');
    csvwrite(strrep(fullFileName,'.tif','_outline.csv'),outline_matrix');
end