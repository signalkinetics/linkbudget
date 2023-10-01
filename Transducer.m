% Transducer.m A class that captures the parameters of an underwater
% acoustic transducer used in acoustic backscatter communications.
classdef Transducer
    properties(Access=public)
        freq % frequency range for transducer properties [Hz]
        Z % Complex electrical impedance of the transducer (R+jX) [Ohm]
        DI % Directivity index of the transducer
        eta % Electromechanical efficiency [1]
        G % Gain [dB]
        TVR % Transmit voltage response [dB]
        role % Role of the transducer in the backscatter system ("tx","rx","node")
        type % String rperesenting transducer type/model/name
        rep % Representation used to define the transducer properties ("Gain", "TVR", "DI & η")
    end
    methods(Access=public)

        % The string in rep defines P1 and P2.
        % Gain: P1=G
        % TVR: P1=TVR
        % DI & η: P1=DI and P2=η
        function tr=Transducer(role,freq,Z,type,rep,P1,P2)
            tr.role=role;
            tr.freq=freq;
            tr.Z=Z;
            tr.type=type;
            tr.rep=rep;
            switch rep
                case "Gain"
                    tr.G=P1;
                    tr.eta=1;
                    tr.DI=tr.G-10*log10(tr.eta);
                    tr.TVR=tr.G+170.8+10*log10(real(tr.Z)./abs(tr.Z).^2);
                case "TVR"
                    tr.TVR=P1;
                    tr.eta=1;
                    tr.G=tr.TVR-10*log10(real(tr.Z)./abs(tr.Z).^2)-170.8;
                    tr.DI=tr.G;
                    tr.DI=tr.G-10*log10(tr.eta);
                case "DI & η"
                    tr.DI=P1;
                    tr.eta=P2;
                    tr.G=tr.DI+10*log10(tr.eta);
                    tr.TVR=tr.G+170.8+10*log10(real(tr.Z)./abs(tr.Z).^2);
            end

            if length(tr.G)==1
                tr.G=tr.G.*ones(size(tr.freq));
            end
        end
        
        
    end
end