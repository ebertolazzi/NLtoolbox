clc;
clear all;
close all;

global nl;

nl = NLtest();

res = nl.listall();
for k=1:length(res)
  fprintf('test N.%03d [%s]\n',k,res{k});
end
