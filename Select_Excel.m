function [A, Text, Row, PipValue, Digit, Dist_HC, Stop_HD] = Select_Excel(k)
%----------------------------------  M5  ----------------------------------
if k == 1
    [A, Text, Row] = xlsread('AUDCHF_M5.csv');
    PipValue = 11;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 2
    [A, Text, Row] = xlsread('AUDJPY_M5.csv');
    PipValue = 6;
    Digit = 3;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 3
    [A, Text, Row] = xlsread('CADJPY_M5.csv');
    PipValue = 6;
    Digit = 3;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 4
    [A, Text, Row] = xlsread('EURNZD_M5.csv');
    PipValue = 6;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 5
    [A, Text, Row] = xlsread('GBPCAD_M5.csv');
    PipValue = 7.33;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 6
    [A, Text, Row] = xlsread('GBPCHF_M5.csv');
    PipValue = 11;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 7
    [A, Text, Row] = xlsread('NZDUSD_M5.csv');
    PipValue = 10;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 8
    [A, Text, Row] = xlsread('XAUUSD_M5.csv');
    PipValue = 10;
    Digit = 2;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
    %---------------
elseif k == 23
    [A, Text, Row] = xlsread('AUDCAD_M5.csv');
    PipValue = 7.33;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 24
    [A, Text, Row] = xlsread('AUDNZD_M5.csv');
    PipValue = 6;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 25
    [A, Text, Row] = xlsread('AUDUSD_M5.csv');
    PipValue = 10;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 26
    [A, Text, Row] = xlsread('CADCHF_M5.csv');
    PipValue = 11;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 27
    [A, Text, Row] = xlsread('CHFJPY_M5.csv');
    PipValue = 6;
    Digit = 3;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 28
    [A, Text, Row] = xlsread('EURAUD_M5.csv');
    PipValue = 6.74;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 29
    [A, Text, Row] = xlsread('EURCAD_M5.csv');
    PipValue = 7.33;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 30
    [A, Text, Row] = xlsread('EURCHF_M5.csv');
    PipValue = 11;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 31
    [A, Text, Row] = xlsread('EURGBP_M5.csv');
    PipValue = 13.38;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 32
    [A, Text, Row] = xlsread('EURJPY_M5.csv');
    PipValue = 6;
    Digit = 3;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 33
    [A, Text, Row] = xlsread('EURUSD_M5.csv');
    PipValue = 10;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 34
    [A, Text, Row] = xlsread('GBPAUD_M5.csv');
    PipValue = 6.74;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 35
    [A, Text, Row] = xlsread('GBPJPY_M5.csv');
    PipValue = 6;
    Digit = 3;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 36
    [A, Text, Row] = xlsread('GBPUSD_M5.csv');
    PipValue = 10;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 37
    [A, Text, Row] = xlsread('NZDCAD_M5.csv');
    PipValue = 7.33;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 38
    [A, Text, Row] = xlsread('NZDCHF_M5.csv');
    PipValue = 11;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 39
    [A, Text, Row] = xlsread('NZDJPY_M5.csv');
    PipValue = 6;
    Digit = 3;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 40
    [A, Text, Row] = xlsread('USDCAD_M5.csv');
    PipValue = 7.33;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 41
    [A, Text, Row] = xlsread('USDCHF_M5.csv');
    PipValue = 11;
    Digit = 5;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 42
    [A, Text, Row] = xlsread('USDJPY_M5.csv');
    PipValue = 6;
    Digit = 3;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
elseif k == 43
    [A, Text, Row] = xlsread('XAGUSD_M5.csv');
    PipValue = 10;
    Digit = 3;
    Dist_HC = 10 * 10^(-(Digit-1));
    Stop_HD = 14 * 10^(-(Digit-1));
%----------------------------------  M1  ----------------------------------
elseif k == 9
    [A, Text, Row] = xlsread('AUDNZD_M1.csv');
    PipValue = 6;
    Digit = 5;
    Dist_HC = 5 * 10^(-(Digit-1));
    Stop_HD = 9 * 10^(-(Digit-1));
elseif k == 10
    [A, Text, Row] = xlsread('AUDUSD_M1.csv');
    PipValue = 10;
    Digit = 5;
    Dist_HC = 5 * 10^(-(Digit-1));
    Stop_HD = 9 * 10^(-(Digit-1));
elseif k == 11
    [A, Text, Row] = xlsread('EURAUD_M1.csv');
    PipValue = 6.74;
    Digit = 5;
    Dist_HC = 5 * 10^(-(Digit-1));
    Stop_HD = 9 * 10^(-(Digit-1));
%----------------------------------  M15  ---------------------------------
elseif k == 12
    [A, Text, Row] = xlsread('AUDCAD_M15.csv');
    PipValue = 7.33;
    Digit = 5;
    Dist_HC = 20 * 10^(-(Digit-1));
elseif k == 13
    [A, Text, Row] = xlsread('AUDCHF_M15.csv');
    PipValue = 11;
    Digit = 5;
    Dist_HC = 20 * 10^(-(Digit-1));
elseif k == 14
    [A, Text, Row] = xlsread('EURCAD_M15.csv');
    PipValue = 7.33;
    Digit = 5;
    Dist_HC = 20 * 10^(-(Digit-1));
elseif k == 15
    [A, Text, Row] = xlsread('EURJPY_M15.csv');
    PipValue = 6;
    Digit = 3;
    Dist_HC = 20 * 10^(-(Digit-1));
elseif k == 16
    [A, Text, Row] = xlsread('EURNZD_M15.csv');
    PipValue = 6;
    Digit = 5;
    Dist_HC = 20 * 10^(-(Digit-1));
elseif k == 17
    [A, Text, Row] = xlsread('GBPJPY_M15.csv');
    PipValue = 6;
    Digit = 3;
    Dist_HC = 20 * 10^(-(Digit-1));
elseif k == 18
    [A, Text, Row] = xlsread('GBPUSD_M15.csv');
    PipValue = 10;
    Digit = 5;
    Dist_HC = 20 * 10^(-(Digit-1));
%----------------------------------  H1  ----------------------------------
elseif k == 19
    [A, Text, Row] = xlsread('AUDNZD_H1.csv');
    PipValue = 6;
    Digit = 5;
    Dist_HC = 30 * 10^(-(Digit-1));
elseif k == 20
    [A, Text, Row] = xlsread('CADCHF_H1.csv');
    PipValue = 11;
    Digit = 5;
    Dist_HC = 30 * 10^(-(Digit-1));
elseif k == 21
    [A, Text, Row] = xlsread('USDCHF_H1.csv');
    PipValue = 11;
    Digit = 5;
    Dist_HC = 30 * 10^(-(Digit-1));
elseif k == 22
    [A, Text, Row] = xlsread('USDJPY_H1.csv');
    PipValue = 6;
    Digit = 3;
    Dist_HC = 30 * 10^(-(Digit-1));
end
end

