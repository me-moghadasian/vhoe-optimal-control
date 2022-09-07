clc

disp('.--------- Vectorised High Order Expansions ---------.')
disp(':          ver. 0.0 (base)                           :')
disp(':....................................................:')
disp('')

disp('.----------------------------------------------------.')
disp(': Vectorised High Order Expansions (VHOE) method is  :')
disp(':    developed by Dr. Mehdi Moghadasian in 2019.     :')
disp(':....................................................:')
disp('')
disp('.----------------------------------------------------.')
% Adding lib path
disp([': Adding path ',pwd,'\lib\']);
addpath([pwd,'\lib\'],'-begin');
disp(': Setup finished.                                    :');
disp(':....................................................:')

% Save paths
savepath;