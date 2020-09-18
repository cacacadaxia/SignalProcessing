

w = f1()
v = f2()
w = f1()
w = f1()

v = f2()
function out = f1()
persistent y;
if isempty(y)
    y = 0;
end
y = y+1;
out = y;
end
function out = f2()
persistent y;
if isempty(y)
    y = 0;
end
y = y + 1;
out = y;
end

