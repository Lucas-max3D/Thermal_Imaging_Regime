% generate ESL for single parameter prints for the thermal imaging study

% first we create the surrounding which is gonna be 16mm by 6.5mm
P_contour=100;
V_contour = 350;
Points= [-3.25,-8,0.05; -3.25, 8,0.05 ; 3.25,8,0.05;3.25,-8,0.05;-3.25, -8,0.05];
Distance_Lines = 0.2;
Distance_Groups = 0.4;
numLinesPerGroup = 5;
numGroupsPerCube = 12;
PPC = importfile('Experimental_parameters.csv');
time_between_lines = 0.03; %delay between two neighboring lines
Start_line_y = -7.5;
length_Line = 6;
layerThickness = 0.05;
z0 = 0.05+layerThickness;

P = PPC(1:12,3);
V = PPC(1:12,2);
ESL = generate_contour(Points,P_contour,V_contour);

%%
ESL = [ESL; generate_group_center_mark(Start_line_y, numLinesPerGroup,numGroupsPerCube, Distance_Lines, Distance_Groups,Points, P_contour,V_contour,ESL(end,3))];
ESL1 = [ESL;generate_all_groups_for_cube(P,V,Distance_Lines,numLinesPerGroup,Distance_Groups,numGroupsPerCube, time_between_lines,Start_line_y,PPC,length_Line,z0,layerThickness)];

P = PPC(13:24,3);
V = PPC(13:24,2);
ESL2 = [ESL;generate_all_groups_for_cube(P,V,Distance_Lines,numLinesPerGroup,Distance_Groups,numGroupsPerCube, time_between_lines,Start_line_y,PPC,length_Line,z0,layerThickness)];

P = PPC(25:36,3);
V = PPC(25:36,2);
ESL3 = [ESL;generate_all_groups_for_cube(P,V,Distance_Lines,numLinesPerGroup,Distance_Groups,numGroupsPerCube, time_between_lines,Start_line_y,PPC,length_Line,z0,layerThickness)];

P = PPC(37:48,3);
V = PPC(37:48,2);
ESL4 = [ESL;generate_all_groups_for_cube(P,V,Distance_Lines,numLinesPerGroup,Distance_Groups,numGroupsPerCube, time_between_lines,Start_line_y,PPC,length_Line,z0,layerThickness)];


PTest = [260 280 320 340];
Vtest = [500 800 1100 1400];
ESLTest = [ESL;generate_all_groups_for_cube(PTest,Vtest,Distance_Lines,numLinesPerGroup,Distance_Groups,numGroupsPerCube, time_between_lines,Start_line_y,PPC,length_Line,z0,layerThickness)];

plotESL(ESLTest)


writeToXML(ESL1,"cube1")
writeToXML(ESL2,"cube2")
writeToXML(ESL3,"cube3")
writeToXML(ESL4,"cube4")
writeToXML(ESLTest,"Test_cube")

%%

function plotESL(ESL)
figure()
for i=1:length(ESL)-1
    if ESL(i,4)>0
        plot(ESL(i:i+1,1),ESL(i:i+1,2),'r')
        hold on
    else
        plot(ESL(i:i+1,1),ESL(i:i+1,2),'b')
        hold on
    end
end
hold off
       

end


function ESL = generate_all_groups_for_cube(P,V,Distance_Lines,numLinesPerGroup,Distance_Groups,numGroupsPerCube, time_between_lines,Start_line_y,PPC,length_Line,z0,layerThickness)
    ESL=[];
    for i=1:length(P)
        ESL = [ESL; generate_Line_Groups(P(i),V(i),Distance_Lines,numLinesPerGroup,Distance_Groups,numGroupsPerCube, time_between_lines,Start_line_y+(i-1)*((numLinesPerGroup-1)*Distance_Lines)+(i-1)*Distance_Groups,PPC,length_Line,z0+(i-1)*layerThickness)];
    end
end



function ESL = generate_Line_Groups(P,V,Distance_Lines,numLinesPerGroup,Distance_Groups,numGroupsPerCube, time_between_lines,Y_start,PPC,length_Line,z)
    k=1;
    for ii=1:numLinesPerGroup
        %this is to account for any potential delays of the laser and the
        %scanhead
        ESL(k,:) = [(-length_Line/2)-(time_between_lines/2)*V,Y_start+ (ii-1)*Distance_Lines,z,0,V];
        k=k+1;
        % This is actual scanning
        ESL(k,:)= [-length_Line/2, Y_start+ (ii-1)*Distance_Lines,z,P,V];
        k=k+1;
        ESL(k,:)= [+length_Line/2, Y_start+ (ii-1)*Distance_Lines,z,0,V];
        k=k+1;
        %once again to compensate for delays and allow enough time for the
        %cooling
        ESL(k,:) = [(length_Line/2)+(time_between_lines/2)*V,Y_start+ (ii-1)*Distance_Lines,z,0,V];
        k=k+1;
    end
end

% function to generate a bounding rectangle for the lines
function ESL = generate_contour(Points, P,V)
    for i=1:length(Points)
        ESL(i,:) = [Points(i,:),P,V];
    end
    ESL(end,4)=0;
end

%This function generates a center mark for each group of lines, this should
%be used to focus the camera between groups
function ESL = generate_group_center_mark(Start_line_y, numLinesPerGroup,numGroupsPerCube, Distance_Lines, Distance_Groups,Points, P,V,z)
    k=1;
    for ii=1:numGroupsPerCube
        ESL(k,:) = [min(Points(:,1)),Start_line_y + ((numLinesPerGroup-1)*Distance_Lines)/2+(ii-1)*(numLinesPerGroup-1)*Distance_Lines+(ii-1)*Distance_Groups,z, P,V];
        k=k+1;
        ESL(k,:) = [min(Points(:,1))-0.4,Start_line_y+ ((numLinesPerGroup-1)*Distance_Lines)/2+(ii-1)*(numLinesPerGroup-1)*Distance_Lines+(ii-1)*Distance_Groups,z, 0,V];
        k=k+1;
        ESL(k,:) = [max(Points(:,1)),Start_line_y+ ((numLinesPerGroup-1)*Distance_Lines)/2+(ii-1)*(numLinesPerGroup-1)*Distance_Lines+(ii-1)*Distance_Groups,z, P,V];
        k=k+1;
        ESL(k,:) = [max(Points(:,1))+0.4,Start_line_y+ ((numLinesPerGroup-1)*Distance_Lines)/2+(ii-1)*(numLinesPerGroup-1)*Distance_Lines+(ii-1)*Distance_Groups,z, 0,V];
        k=k+1;
    end
end

function Experimentalparameters = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  EXPERIMENTALPARAMETERS = IMPORTFILE(FILENAME) reads data from text
%  file FILENAME for the default selection.  Returns the numeric data.
%
%  EXPERIMENTALPARAMETERS = IMPORTFILE(FILE, DATALINES) reads data for
%  the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  Experimentalparameters = importfile("C:\Users\schlenge\OneDrive - epfl.ch\PhD\Thermal_imaging_regime\Experimental_parameters.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 17-Oct-2024 16:08:33

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "Velocity", "Power", "Order"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
Experimentalparameters = readtable(filename, opts);

%% Convert to output type
Experimentalparameters = table2array(Experimentalparameters);
end

function writeToXML(ESL1,name)
LayerIdx= find(diff(ESL1(:,3))>0);
LayerIdx= [0;LayerIdx;length(ESL1)];
ESL_Global={};
for i=1:length(LayerIdx)-1
    ESL_Global{i}= ESL1((LayerIdx(i)+1):LayerIdx(i+1),:);
end

Thickness = 0.05; %layer thickness in mm
mkdir(name);
tV1=0;
tV2 =0;
tL1=0;
tL2 =0;
SyncDelay=0;
PartNo =0;
jump_speed = 4000;
% lets loop over the layers
for i=1:size(ESL_Global,2)
    number = convertCharsToStrings(sprintf( '%03d', i )) ;
    fileID = fopen(fullfile(name,number+'.xml'),'w');
    % first I guess we need to find all velocities.
    Velocities = unique(ESL_Global{1,i}(:,5));
    
    % write the first build and the layer thickness
    fprintf(fileID,'%s\n','<Build>');
    fprintf(fileID,'\t%s','<Thickness>');
    fprintf(fileID,'%.5f',Thickness);
    fprintf(fileID,'%s\n','</Thickness>');

    % start the velocitylist
    fprintf(fileID,'\t%s\n','<VelocityProfileList>')

    for j=1:length(Velocities)
        fprintf(fileID,'\t\t%s\n','<VelocityProfile>')
        
        fprintf(fileID,'\t\t\t%s','<ID>');
        fprintf(fileID,'%03d',j);
        fprintf(fileID,'%s\n','</ID>');
        
        fprintf(fileID,'\t\t\t%s','<Velocity>');
        fprintf(fileID,'%03d',Velocities(j));
        fprintf(fileID,'%s\n','</Velocity>');
        if Velocities(j)==jump_speed
            mode='Jump';
        else
            mode='Hatching-Infill';
        end
        fprintf(fileID,'\t\t\t%s','<Mode>');
        fprintf(fileID,'%s',mode);
        fprintf(fileID,'%s\n','</Mode>');
        
        %Delays
        fprintf(fileID,'\t\t\t%s','<tV1>');
        fprintf(fileID,'%d',tV1);
        fprintf(fileID,'%s\n','</tV1>');

        fprintf(fileID,'\t\t\t%s','<tV2>');
        fprintf(fileID,'%d',tV2);
        fprintf(fileID,'%s\n','</tV2>');

        fprintf(fileID,'\t\t\t%s','<tL1>');
        fprintf(fileID,'%d',tL1);
        fprintf(fileID,'%s\n','</tL1>');

        fprintf(fileID,'\t\t\t%s','<tL2>');
        fprintf(fileID,'%d',tL2);
        fprintf(fileID,'%s\n','</tL2>');
        fprintf(fileID,'\t\t%s\n','</VelocityProfile>')
    end
    fprintf(fileID,'\t%s\n','</VelocityProfileList>')

        % now over to the actuall scanning path 
        fprintf(fileID,'\t%s\n','<Trajectory>')
        fprintf(fileID,'\t\t%s\n','<TravelerID>0</TravelerID>')
        
        fprintf(fileID,'\t\t%s','<SyncDelay>');
        fprintf(fileID,'%d',SyncDelay);
        fprintf(fileID,'%s\n','</SyncDelay>');

        fprintf(fileID,'\t\t%s\n','<Path>')
        fprintf(fileID,'\t\t\t%s\n','<Type>Hatch</Type>')

        fprintf(fileID,'\t\t\t%s','<Tag>');
        fprintf(fileID,'%s','Part ');
        fprintf(fileID,'%d;',PartNo);
        fprintf(fileID,'%s','Layer ');
        fprintf(fileID,'%d;',i);
        fprintf(fileID,'%s','Power ');
        fprintf(fileID,'%d;',max(ESL_Global{1,i}(:,4)));
        fprintf(fileID,'%s','Velocity ');
        fprintf(fileID,'%d;',max(ESL_Global{1,i}(:,5)));
        fprintf(fileID,'%s;','Hatch');
        fprintf(fileID,'%s\n','</Tag>');

        %now the start point
        numSegments = length(ESL_Global{1,i}(:,5))-1;
        
        fprintf(fileID,'\t\t\t%s','<NumSegments>');
        fprintf(fileID,'%d',numSegments);
        fprintf(fileID,'%s\n','</NumSegments>');

        fprintf(fileID,'\t\t\t%s\n','<Start>');
        
        fprintf(fileID,'\t\t\t\t%s','<X>');
        fprintf(fileID,'%0.6f',ESL_Global{1,i}(1,1));
        fprintf(fileID,'%s\n','</X>');
        fprintf(fileID,'\t\t\t\t%s','<Y>');
        fprintf(fileID,'%0.6f',ESL_Global{1,i}(1,2));
        fprintf(fileID,'%s\n','</Y>');
        fprintf(fileID,'\t\t\t%s\n','</Start>');

        % followed by the segments
        for j=2:length(ESL_Global{1,i}(:,1))
            
            X = ESL_Global{1,i}(j,1);
            Y = ESL_Global{1,i}(j,2);
            Power = ESL_Global{1,i}(j-1,4);
            VelocityID =find(Velocities==ESL_Global{1,i}(j-1,5));
            
            fprintf(fileID,'\t\t\t%s\n','<Segment>');
            fprintf(fileID,'\t\t\t\t%s','<SegmentID>');
            fprintf(fileID,'%d',j-1);
            fprintf(fileID,'%s\n','</SegmentID>');
            
            fprintf(fileID,'\t\t\t\t%s','<Power>');
            fprintf(fileID,'%d',Power);
            fprintf(fileID,'%s\n','</Power>');

            fprintf(fileID,'\t\t\t\t%s','<idxVelocityProfile>');
            fprintf(fileID,'%03d',VelocityID);
            fprintf(fileID,'%s\n','</idxVelocityProfile>');

            fprintf(fileID,'\t\t\t\t%s\n','<End>');

            fprintf(fileID,'\t\t\t\t\t%s','<X>');
            fprintf(fileID,'%0.6f',X);
            fprintf(fileID,'%s\n','</X>');
            fprintf(fileID,'\t\t\t\t\t%s','<Y>');
            fprintf(fileID,'%0.6f',Y);
            fprintf(fileID,'%s\n','</Y>');

            fprintf(fileID,'\t\t\t\t%s\n','</End>');
            fprintf(fileID,'\t\t\t%s\n','</Segment>');




        end
        fprintf(fileID,'\t\t%s\n','</Path>');
        fprintf(fileID,'\t%s\n','</Trajectory>');
        fprintf(fileID,'%s','</Build>');
        fclose(fileID)
end

end

