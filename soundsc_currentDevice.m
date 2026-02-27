% This function can be called instead of 'soundsc' if you want to change the output device during a
% Matlab session (without having to restart Matlab)
function soundsc_currentDevice(y, fs)
% Get the current macOS default output device ID
deviceID = get_default_output_deviceID();  % 0 = output, returns default device ID

% Normalise like soundsc does
y = y / max(abs(y(:)));

p = audioplayer(y, fs, 16, deviceID);
play(p);
waitfor(p, 'Running', 'off');

    function ID = get_default_output_deviceID()
        % Ask macOS for the default output device name
        [~, devname] = system('/opt/homebrew/bin/SwitchAudioSource -c');

        % Match against MATLAB's device list
        qq          = audiodevinfo;

        [~,idx] = find(strcmp({qq.output.Name}, char([devname(1:end-1) ' (Core Audio)'])));
        ID = qq.output(idx).ID;
    end
end