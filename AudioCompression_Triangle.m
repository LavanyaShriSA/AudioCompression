
clear;
close all;

% Read the audio data from file
[audio_data fs] = audioread("cleanspeech.wav");

% Create an output file to store SNR values
file_id = fopen('project_2_triangle.csv','w');
fprintf(file_id, "N,n,Method 1 overall SNR,Method 1 segmental SNR,Method 2 overall SNR,Method 2 segmental SNR\n");

N_values = [64, 128, 256];
lines = ["r", "b", "g"];

ax_seg_1 = subplot(1, 1, 1, 'Parent', figure());
ax_seg_2 = subplot(1, 1, 1, 'Parent', figure());
ax_ovr_1 = subplot(1, 1, 1, 'Parent', figure());
ax_ovr_2 = subplot(1, 1, 1, 'Parent', figure());

% Repeat experiment for different N values
for i=1:length(N_values)

    N = N_values(i);

    % Create a triangular window to multiple with each frame
    triang_window = triang(N);

    segmental_snrs_method_1 = [];
    segmental_snrs_method_2 = [];
    overall_snrs_method_1 = [];
    overall_snrs_method_2 = [];
    x_values = [];

    data_length = length(audio_data);

    % vary n from 1 to 1+(N/2)
    for n=1:(1+(N/2))

        frames_reconstruct_1 = [];
        frames_reconstruct_2 = [];

        % Divide data into frame of length N with overlap of N/2
        for start=1:N/2:(data_length-(N/2))

            stop = start+N-1;

            % When last frame length is less than N pad zeros to length N
            if stop > data_length
                frame_data = audio_data(start:data_length, :);
                frame_data = cat(1, frame_data, zeros(stop-data_length, 1));
            else
                frame_data = audio_data(start:stop, :);
            end

            % Calulate the FFT after multiplying the data with triangular window
            frame_fft = fft(frame_data.*triang_window, N);

            % Method 1 - dominant n values

            % Sort the values and use only n dominant (large) values
            % Use only 1+(N/2) values due to symmetry of FFT
            [sorted_values,sorted_index] = sort(abs(frame_fft(1:1+(N/2))),'descend');
            dominant_index = sorted_index(1:n);

            % Use zeros in places of non dominant values
            dominant_fft = zeros(N,1);

            % Fill the dominant values
            for idx=1:length(dominant_index)
                idx=dominant_index(idx);

                dominant_fft(idx) = frame_fft(idx);

                % Fill the other side with domianant values due to symmetry
                % Ex [a,0,b,c,0,0,c*,b*,0,a*]
                if N-idx+2<=N
                    dominant_fft(N-idx+2) = frame_fft(N-idx+2);
                end
            end

            % Reconstruct using inverse fft
            x_reconstruct1 = ifft(dominant_fft);

            if start==1
                frames_reconstruct_1 = x_reconstruct1;
            else

                % Concatenate the reconstructed sequence with previous sequence along with overlap by padding with zeros
                % previous sequence  [x x x x x x 0 0 0]
                % Current sequence   [0 0 0 y y y y y y]

                frames_reconstruct_1 = vertcat(frames_reconstruct_1, zeros(N/2,1));
                x_reconstruct1 = vertcat(zeros(start-1,1), x_reconstruct1);
                frames_reconstruct_1 = frames_reconstruct_1 + x_reconstruct1;
            end


            % Method 2 - First n values

            % Use only first n values of FFT and fill others with 0
            % Due to symmetry of FFT, complex conjucates of first n values from other side is also included
            first_fft = frame_fft;
            for idx=n+1:N-n+1
                first_fft(idx) = 0;
            end

            % Reconstruct using inverse fft
            x_reconstruct2 = ifft(first_fft);

            if start==1
                frames_reconstruct_2 = x_reconstruct2;
            else
                % Concatenate the reconstructed sequence with previous sequence along with overlap by padding with zeros
                % previous sequence  [x x x x x x 0 0 0]
                % Current sequence   [0 0 0 y y y y y y]

                frames_reconstruct_2 = vertcat(frames_reconstruct_2, zeros(N/2,1));
                x_reconstruct2 = vertcat(zeros(start-1,1), x_reconstruct2);
                frames_reconstruct_2 = frames_reconstruct_2 + x_reconstruct2;
            end



        end

        % Calculate segmental and overall SNR for method 1
        [segmental_snr_1, overall_snr_1] = calculate_snr(audio_data, frames_reconstruct_1, N);

        % Calculate segmental and overall SNR for method 2
        [segmental_snr_2, overall_snr_2] = calculate_snr(audio_data, frames_reconstruct_2, N);

        % Store the calculated SNR values
        segmental_snrs_method_1 = [segmental_snrs_method_1 segmental_snr_1];
        segmental_snrs_method_2 = [segmental_snrs_method_2 segmental_snr_2];
        overall_snrs_method_1 = [overall_snrs_method_1 overall_snr_1];
        overall_snrs_method_2 = [overall_snrs_method_2 overall_snr_2];

        % Store the n/N value for plotting
        x_values = [x_values n/N];

        % Write the SNR values in output file
        fprintf(file_id, "%d,%d,%f,%f,%f,%f\n", N, n, overall_snr_1, segmental_snr_1, overall_snr_2, segmental_snr_2);

    end

    % Plot the segmental SNR for method 1
    plot(ax_seg_1, x_values, segmental_snrs_method_1, lines(i));
    hold(ax_seg_1, "on");

    % Plot the segmental SNR for method 2
    plot(ax_seg_2, x_values, segmental_snrs_method_2, lines(i));
    hold(ax_seg_2, "on");

    % Plot the overall SNR for method 1
    plot(ax_ovr_1, x_values, overall_snrs_method_1, lines(i));
    hold(ax_ovr_1, "on");

    % Plot the overall SNR for method 2
    plot(ax_ovr_2, x_values, overall_snrs_method_2, lines(i));
    hold(ax_ovr_2, "on");

end

% Close the output file
fclose(file_id);

% Add labels, title and legend for the plots
legend(ax_seg_1, ["N=64", "N=128", "N=256"], "Location", "north");
title(ax_seg_1, "Segmental SNR - Method 1");
xlabel(ax_seg_1, "n/N ratio");
ylabel(ax_seg_1, "Segmental SNR (dB)");
grid(ax_seg_1, "on");

legend(ax_seg_2, ["N=64", "N=128", "N=256"], "Location", "north");
title(ax_seg_2, "Segmental SNR - Method 2");
xlabel(ax_seg_2, "n/N ratio");
ylabel(ax_seg_2, "Segmental SNR (dB)");
grid(ax_seg_2, "on");

legend(ax_ovr_1, ["N=64", "N=128", "N=256"], "Location", "north");
title(ax_ovr_1, "Overall SNR - Method 1");
xlabel(ax_ovr_1, "n/N ratio");
ylabel(ax_ovr_1, "Overall SNR (dB)");
grid(ax_ovr_1, "on");

legend(ax_ovr_2, ["N=64", "N=128", "N=256"], "Location", "north");
title(ax_ovr_2, "Overall SNR - Method 2");
xlabel(ax_ovr_2, "n/N ratio");
ylabel(ax_ovr_2, "Overall SNR (dB)");
grid(ax_ovr_2, "on");

function [seg_snr, overall_snr] = calculate_snr(x_orig, x_reconstructed, N)

    snr_num = 0;
    snr_den = 0;
    seg_snr = 0;
    shape = size(x_orig);
    L = 0;

    for frame_no = 1:N:shape(1)
        x_frame = x_orig(frame_no:frame_no+N-1, :);
        error = x_orig(frame_no:frame_no+N-1, :) - x_reconstructed(frame_no:frame_no+N-1, :);

        num = (transpose(x_frame)*x_frame);
        den = (transpose(error)*error);

        % Calculate the SNR for each frame
        seg_snr = seg_snr + log10( num / den );

        % Update the numerate and denominator for each frame
        snr_num = snr_num + num;
        snr_den = snr_den + den;

        L = L + 1;

    end

    % Caclulate the overall SNR
    overall_snr = 10*log10(snr_num/snr_den);

    % Find the average SNR of all frames
    seg_snr = (10/L) * seg_snr;

end




