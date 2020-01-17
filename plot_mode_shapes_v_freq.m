freqs = logspace(0, 6,120);
num_points = 100;
beam_points = linspace(0,.15,num_points);
mode_shapes = zeros(length(freqs), num_points);

driving_force = 1;

tic
for index = 1:length(freqs)
    mode_shapes(index,:) = mode_shape_from_freq(freqs(index),...
                                  num_points, driving_force);
    disp(index);
end
toc

abs_mode_shapes = abs(mode_shapes);
db_mode_shapes = 20*log10(abs_mode_shapes);

abs_accel_shapes = (((2*pi*freqs).^2)  .*(abs_mode_shapes.')).';
db_accel_shapes  = 20*log10(abs_accel_shapes);

figure();
surf(beam_points, freqs, db_mode_shapes);
set(gca, 'YScale', 'log')
xlabel("beam distance (m)")
ylabel("freq (Hz)")
zlabel("response (~db)")
title("Beam deflection response mode vs driving frequency")

figure();
surf(beam_points, freqs, db_accel_shapes);
set(gca, 'YScale', 'log')
xlabel("beam distance (m)")
ylabel("freq (Hz)")
zlabel("response (~db)")
title("Beam acceleration response mode vs driving frequency")

figure();
semilogx(freqs, db_mode_shapes(:,18))
title("Displacement Frequency response at 1 inch");
xlabel("Freq (Hz)")
ylabel("Deflection Response (dB) (m/N)")

figure();
semilogx(freqs, db_accel_shapes(:,18))
title("acceleration Frequency response at 1 inch");
xlabel("Freq (Hz)")
ylabel("Accel Response (dB) ((m/s^2)/N)")
