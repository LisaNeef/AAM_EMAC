
%And my data is plotted in a loop so I have

for i = 1:6
h(i) = subplot(3,2,i)
etc.
end

%Then this after the loop (so it plots everything and then adjusts it). Anyway, that is why I have h(1) instead of h1.


set(gcf,'Position',[100 65 572 719])

x0 = 0.09;
y0 = 0.99;
w = 0.45;
w2 = w-0.01;
ht = y0/3;

ht2 = ht;
y3 = y0-ht+0;
set(h(1),'Position',[x0 y3 w2 ht2])
set(h(2),'Position',[x0+w y3 w2 ht2])

ht2 = ht;
y3 = y0-2*ht+0.05;
set(h(3),'Position',[x0 y3 w2 ht2])
set(h(4),'Position',[x0+w y3 w2 ht2])

ht2 = ht;
y3 = y0-3*ht+0.16;
set(h(5),'Position',[x0 y3 w2 ht2])
set(h(6),'Position',[x0+w y3 w2 ht2])
