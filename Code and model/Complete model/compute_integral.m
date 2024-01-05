function I = compute_integral(t, y)

I = 0; 

for i = 1 : length(t) - 1

    I = I + (y(i + 1) + y(i) ) * (t(i + 1) - t(i)) * 0.5;

end

end