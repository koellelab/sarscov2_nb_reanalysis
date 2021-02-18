function void = main_readInData(void)

clear all; close all; clc;

for i = 1:39
    
    switch i
        case 1
            values = readmatrix('CoV_003_CoV_143_0.01.csv');
            data(i).donor = 003; data(i).recipient = 143;
        case 2
            values = readmatrix('CoV_003_CoV_146_0.01.csv');
            data(i).donor = 003; data(i).recipient = 146;
        case 3
            values = readmatrix('CoV_146_CoV_147_0.01.csv');
            data(i).donor = 146; data(i).recipient = 147;
        case 4
            values = readmatrix('CoV_146_CoV_172_0.01.csv');
            data(i).donor = 146; data(i).recipient = 172;
        case 5
            values = readmatrix('CoV_146_CoV_157_0.01.csv');
            data(i).donor = 146; data(i).recipient = 157;
        case 6
            values = readmatrix('CoV_146_CoV_159_0.01.csv');
            data(i).donor = 146; data(i).recipient = 159;
        case 7
            values = readmatrix('CoV_146_CoV_150_0.01.csv');
            data(i).donor = 146; data(i).recipient = 150;
        case 8
            values = readmatrix('CoV_146_CoV_1056_0.01.csv');
            data(i).donor = 146; data(i).recipient = 1056;
        case 9
            values = readmatrix('CoV_159_CoV_168_0.01.csv');
            data(i).donor = 159; data(i).recipient = 168;
        case 10
            values = readmatrix('CoV_159_CoV_169_0.01.csv');
            data(i).donor = 159; data(i).recipient = 169;
        case 11
            values = readmatrix('CoV_150_CoV_162_0.01.csv');
            data(i).donor = 150; data(i).recipient = 162;    
        case 12
            values = readmatrix('CoV_162_CoV_161_0.01.csv');
            data(i).donor = 162; data(i).recipient = 161;
        case 13
            values = readmatrix('CoV_1056_CoV_166_0.01.csv');
            data(i).donor = 1056; data(i).recipient = 166;
        case 14
            values = readmatrix('CoV_1056_CoV_176_0.01.csv');
            data(i).donor = 1056; data(i).recipient = 176;
        case 15
            values = readmatrix('CoV_1056_CoV_171_0.01.csv');
            data(i).donor = 1056; data(i).recipient = 171;
        case 16
            values = readmatrix('CoV_1056_CoV_1057_0.01.csv');
            data(i).donor = 1056; data(i).recipient = 1057;
        case 17
            values = readmatrix('CoV_1056_CoV_1067_0.01.csv');
            data(i).donor = 1056; data(i).recipient = 1067;
        case 18
            values = readmatrix('CoV_171_CoV_180_0.01.csv');
            data(i).donor = 171; data(i).recipient = 180;
        case 19
            values = readmatrix('CoV_180_CoV_197_0.01.csv');
            data(i).donor = 180; data(i).recipient = 197;
        case 20
            values = readmatrix('CoV_1067_CoV_183_0.01.csv');
            data(i).donor = 1067; data(i).recipient = 183;
        case 21
            values = readmatrix('CoV_1057_CoV_175_0.01.csv');
            data(i).donor = 1057; data(i).recipient = 175;
        case 22
            values = readmatrix('CoV_1057_CoV_177_0.01.csv');
            data(i).donor = 1057; data(i).recipient = 177;
        case 23
            values = readmatrix('CoV_1057_CoV_1058_0.01.csv');
            data(i).donor = 1057; data(i).recipient = 1058;
        case 24
            values = readmatrix('CoV_1057_CoV_187_0.01.csv');
            data(i).donor = 1057; data(i).recipient = 187;
        case 25
            values = readmatrix('CoV_1058_CoV_1059_0.01.csv');
            data(i).donor = 1058; data(i).recipient = 1059;
        case 26
            values = readmatrix('CoV_1058_CoV_1064_0.01.csv');
            data(i).donor = 1058; data(i).recipient = 1064;    
        case 27    
            values = readmatrix('CoV_1058_CoV_1063_0.01.csv');
            data(i).donor = 1058; data(i).recipient = 1063;
        case 28
            values = readmatrix('CoV_1059_CoV_1065_0.01.csv');
            data(i).donor = 1059; data(i).recipient = 1065;
        case 29
            values = readmatrix('CoV_1059_CoV_1062_0.01.csv');
            data(i).donor = 1059; data(i).recipient = 1062;
        case 30
            values = readmatrix('CoV_1062_CoV_218_0.01.csv');
            data(i).donor = 1062; data(i).recipient = 218;
        case 31
            values = readmatrix('CoV_1062_CoV_219_0.01.csv');
            data(i).donor = 1062; data(i).recipient = 219;
        case 32 
            values = readmatrix('CoV_1062_CoV_217_0.01.csv');
            data(i).donor = 1062; data(i).recipient = 217;
        case 33
            values = readmatrix('CoV_217_CoV_256_0.01.csv');
            data(i).donor = 217; data(i).recipient = 256;
        case 34
            values = readmatrix('CoV_187_CoV_1068_0.01.csv');
            data(i).donor = 187; data(i).recipient = 1068;
        case 35
            values = readmatrix('CoV_194_CoV_193_0.01.csv');
            data(i).donor = 194; data(i).recipient = 193;
        case 36
            values = readmatrix('CoV_194_CoV_195_0.01.csv');
            data(i).donor = 194; data(i).recipient = 195;
        case 37
            values = readmatrix('CoV_194_CoV_196_0.01.csv');
            data(i).donor = 194; data(i).recipient = 196;
        case 38
            values = readmatrix('CoV_198_CoV_230_0.01.csv');
            data(i).donor = 198; data(i).recipient = 230;
        case 39
            values = readmatrix('CoV_273_CoV_271_0.01.csv');
            data(i).donor = 273; data(i).recipient = 271;      
        otherwise
            error('not a valid transmission pair');
    end
    data(i).donor_iSNVs = values(:,1);
    data(i).recipient_iSNVs = values(:,2);
end

n_TPs = length(data);

sample_list = [];
for i = 1:n_TPs
   sample_list = [sample_list data(i).donor data(i).recipient];
end
unique_sample_list = unique(sample_list);

CT_data.sample_name = unique_sample_list;
CT_data.CT_value = [19.2 24.2 16.8 34.8 24.9 21.4 21.7 21.55 19.92 NaN 24.78 18.86 20.11 28.58 26.87 18.79 28.33 20.44 19.4 34.8 24.2 18.1 22.7 24.5 30.9 22.2 15.59 14.69 16.33 16.6 19.7 31 34.2 19.5 13.9 17.2 16.8 17 28.3 16.1 18 26.7 17.9];

save('Popa_data', 'data', 'n_TPs', 'CT_data');