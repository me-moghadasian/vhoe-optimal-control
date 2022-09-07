clc

disp('.--- Vectorised High Order Expansions Optimal Control ---.')
disp(':    ver. 0.0                                            :')
disp(':........................................................:')
disp('')

disp('.--------------------------------------------------------.')
disp(': Vectorised High Order Expansions (VHOE) method was     :')
disp(':    developed by Dr. Mehdi Moghadasian in 2019.         :')
disp(':........................................................:')
disp('')
disp('.--------------------------------------------------------.')
% Adding lib path
disp([': Adding path ',pwd,'\lib\']);
addpath([pwd,'\lib\'],'-begin');
disp(': Setup finished.                                        :');
disp(':........................................................:')

% Save paths
savepath;