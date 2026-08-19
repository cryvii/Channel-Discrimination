function Ch = draw_one_channel(channel_type, d)
% Draw one Haar-random channel of the requested type.
    if strcmpi(channel_type, 'unitary')
        U  = RandomUnitary(d);
        Ch = ChoiMatrix({U});
    elseif strcmpi(channel_type, 'superoperator') || strcmpi(channel_type, 'channel') || strcmpi(channel_type, 'channels')
        % Corrected terminology mapping:
        % Maps high-level user 'channel' or 'channels' to general Haar channels.
        Ch = RandomSuperoperator(d);
    elseif strcmpi(channel_type, 'pauli')
        if d ~= 2
            error('PauliChannel only supports qubits (d=2).');
        end
        Ch = PauliChannel();  % random theta, random q, no printout suppression
    else
        error('Unknown channel_type: %s', channel_type);
    end
end