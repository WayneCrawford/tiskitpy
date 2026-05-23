% PLOTOUTPUTS Script to plot ouputs of FIR_CORR
%
% W Crawford 3/2013

figure(1); clf
td=   ts_data(load('sgn.bhz')',20);          td.channel(1).metadata.component='BHZ';
td=td+ts_data(load('sgn.bhz.corr.old')',20); td.channel(2).metadata.component='OLD';
td=td+ts_data(load('sgn.bhz.corr')',20);     td.channel(3).metadata.component='NEW';
td=td+ts_data(load('sgn.hhz')',80);          td.channel(4).metadata.component='HHZ';
td=td+ts_data(load('sgn.hhz.corr.old')',80); td.channel(5).metadata.component='OLD';
td=td+ts_data(load('sgn.hhz.corr')',80);     td.channel(6).metadata.component='NEW';

td 

plot(td,'0000/01/01 00:00:17',3);