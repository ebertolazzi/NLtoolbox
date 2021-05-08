clc;
clear all;
close all;

global nl;

nl = NLtest();

N = nl.numTest();
for k=1:N
  fprintf('test N.%03d [%s]\n',k,nl.name(k));
  fprintf('%s\n\n',nl.bibtex(k));
end
