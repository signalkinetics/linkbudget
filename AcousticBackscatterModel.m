classdef AcousticBackscatterModel < handle
    properties (Access=public)
        tx % Tx transducer properties
        rx % Rx transducer properties
        node % Node transducer properties
        channel % Channel properties
        freq=logspace(1,6,1e4) % Default frequency vector
        P_ref=1e-6 % Reference pressure for SPL underwater (1 uPa)
        
    end

    methods (Access = public)
        function model=AcousticBackscatterModel()
            model.loadChannelDefaults();
            model.createSphericalTransducer("tx");
            model.createSphericalTransducer("rx");
            model.createSphericalTransducer("node");
        end

        function loadChannelDefaults(model)
            model.channel.c= 1480;
            model.channel.rho=1000;
            model.channel.k_s=2;
            model.channel.shipping=0.5;
            model.channel.wind=10;
            model.channel.salinity=35;
            model.channel.pH=8;
            model.channel.temp=5;
            model.channel.depth=0;
        end

        function createSphericalTransducer(model,role)
            Z=50+0j;
            rep="Gain";
            type="Spherical";
            model.setTr(Transducer(role,model.freq,Z,type,rep,0));
        end

        function setTr(model,tr)
            switch tr.role
                case "tx"
                    model.tx=tr;
                case "rx"
                    model.rx=tr;
                case "node"
                    model.node=tr;
            end

        end

        function tr=getTr(model,role)
            switch role
                case "tx"
                    tr=model.tx;
                case "rx"
                    tr=model.rx;
                case "node"
                    tr=model.node;
            end
        end

        function loadTransducer(model,role,path,transducerName)
            data=readtable(path+"/"+transducerName+".csv");
            f=data.(1).';
            TVR=data.(2).';
            if contains(upper(data.Properties.VariableNames{3}),["G", "COND"])
                Sig=data.(3)+1j*data.(4);
                Z=(1./Sig).';
            elseif contains(upper(data.Properties.VariableNames{3}),["Z", "IMP"])
                Z=data.(3)*exp(1j*data.(4)/180*pi);
                Z=Z.';
            elseif contains(upper(data.Properties.VariableNames{3}),["SIGMA", "ADM"])
                Sig=data.(3)*exp(1j*data.(4)/180*pi);
                Z=(1./Sig).';
            else
                Z=data.(3).'+1j*data.(4).';
            end
            model.setTr(Transducer(role,f,Z,transducerName,"TVR",TVR));
        end

        function freqRange = getFreqRange(model)
            txmin=min(model.tx.freq);
            txmax=max(model.tx.freq);
            rxmin=min(model.rx.freq);
            rxmax=max(model.rx.freq);
            nodemin=min(model.node.freq);
            nodemax=max(model.node.freq);
            freqRange=[max([txmin rxmin nodemin]) min([txmax rxmax nodemax])];
        end

        function [G_tx,G_rx,G_node,lambda,alpha]=getDerivedParameters(model,freq)
            G_tx=interp1(model.tx.freq,model.tx.G,freq);
            G_rx=interp1(model.rx.freq,model.rx.G,freq);
            G_node=interp1(model.node.freq,model.node.G,freq);
            
            alpha=waterAtten(freq,model.channel.temp,model.channel.salinity,model.channel.depth/1000,model.channel.pH);
            lambda=model.channel.c./freq;

        end

        function SNR=calculateSNR(model,freq,range,Pin,BW,Ze1,Ze2)
            [G_tx,G_rx,G_node,lambda,alpha]=model.getDerivedParameters(freq);
            PL=10*model.channel.k_s*log10(range)+alpha*range;

            if nargin>5
                IML_U=20*log10(real(model.node.Z).*(abs(Ze1 - Ze2)./(abs(Ze1+model.node.Z).*abs(Ze2 + model.node.Z))));
            else
                IML_U=0;
            end
            % Calculate the backscatter level for the uplink
            % channel (subtracting 3 dB to account for energy lost
            % to higher order harmonics due to square wave
            % switching)
            BL=159.8+10*log10(Pin)+G_tx+2*G_node+G_rx+10*log10(lambda.^2/pi) ...
            -2*PL+IML_U - 3;
            NL = noiseDensity(freq/1000,model.channel.shipping,model.channel.wind);

            SNR=BL-(NL+10*log10(BW));
        end

        function P_h=calculateHarvestedPower(model,freq,range,Pin,Ze)
            
            [G_tx,~,G_node,~,alpha]=model.getDerivedParameters(freq);
            PL=10*model.channel.k_s*log10(range)+alpha*range;


            if nargin>4
                IML_D=10*log10((real(Ze).*real(model.node.Z))./(abs(model.node.Z+Ze).^2));
            else
                IML_D=0;
            end
            P_h_dBW = 47.7+10*log10(Pin)+G_tx-PL+G_node- 20*log10(freq) + IML_D; 
            P_h=10.^(P_h_dBW/10);
        end

        function Pin_d=calculateMinPowerInDownlink(model,freq,range,P_th,Ze)
            [G_tx,~,G_node,~,alpha]=model.getDerivedParameters(freq);
            PL=10*model.channel.k_s*log10(range)+alpha*range;

            if nargin>4
                IML_D=10*log10((real(Ze).*real(model.node.Z))./(abs(model.node.Z+Ze).^2));
            else
                IML_D=0;
            end

            Pin_d=10.^((10*log10(P_th)-(170.8+G_tx-PL+G_node...
                    +10*log10(model.channel.c/(model.channel.rho*pi)) ... 
                    + 20*log10(model.P_ref) - 20*log10(freq) + IML_D))/10);                   
        end

        function Pin_u=calculateMinPowerInUplink(model,freq,range,BW,SNR_th,Ze1,Ze2s)
            [G_tx,G_rx,G_node,lambda,alpha]=model.getDerivedParameters(freq);
            PL=10*model.channel.k_s*log10(range)+alpha*range;

            if nargin>5
                IML_U=20*log10(real(model.node.Z).*(abs(Ze1 - Ze2)./(abs(Ze1+model.node.Z).*abs(Ze2 + model.node.Z))));
            else
                IML_U=0;
            end

            NL = noiseDensity(freq/1000,model.channel.shipping,model.channel.wind);

            Pin_u=10.^((SNR_th-(159.8+G_tx+2*G_node+G_rx+10*log10(lambda.^2/pi) ...
                    -2*PL+IML_U - 3-(NL+10*log10(BW))))/10);
        end

        function range_d=calculateMaxRangeDownlink(model,freq,Pin,P_th,Ze)
            [G_tx,~,G_node,~,alpha]=model.getDerivedParameters(freq);
            
            if nargin>4
                IML_D=10*log10((real(Ze).*real(model.node.Z))./(abs(model.node.Z+Ze).^2));
            else
                IML_D=0;
            end

            %Solution for ax+blog10(x)+c=0 is
            %x=(b/a/ln(10))LambertW(a*ln(10)/b*exp(-c*ln10/b))
            a=alpha;
            b=10*model.channel.k_s;
            c=(10*log10(P_th)-(170.8+10*log10(Pin)+G_tx+G_node...
                +10*log10(model.channel.c/(model.channel.rho*pi)) ... 
                + 20*log10(model.P_ref) - 20*log10(freq) + IML_D));
            range_d=(b./a./log(10)).*Lambert_W(a.*log(10)./b.*exp(-c.*log(10)./b));
        end

        function range_u=calculateMaxRangeUplink(model,freq,Pin,SNR_th,BW,Ze1,Ze2)
            %Solution for ax+blog10(x)+c=0 is
            %x=(b/a/ln(10))LambertW(a*ln(10)/b*exp(-c*ln10/b))
            [G_tx,G_rx,G_node,lambda,alpha]=model.getDerivedParameters(freq);
            NL = noiseDensity(freq/1000,model.channel.shipping,model.channel.wind);
            if nargin>5
                IML_U=20*log10(real(model.node.Z).*(abs(Ze1 - Ze2)./(abs(Ze1+model.node.Z).*abs(Ze2 + model.node.Z))));
            else
                IML_U=0;
            end
            
            a=2*alpha;
            b=2*10*model.channel.k_s;
            c=SNR_th-(159.8+10*log10(Pin)+G_tx+2*G_node+G_rx+10*log10(lambda.^2/pi) ...
            +IML_U - 3-(NL+10*log10(BW)));
            range_u=(b./a./log(10)).*Lambert_W(a.*log(10)./b.*exp(-c.*log(10)./b));
        end
    end

end
