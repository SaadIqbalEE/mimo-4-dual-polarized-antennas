clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The QPSK modulation %%%  

tic;        %Timer Started

%Two modulation Schemes QPSK & QPSKr
% QPSK <--- Z1
% QPSKr<--- Z2

l= floor(100000/8)*8;            %No. of Bits @Note: should be divisble by nbits

snrdb=-4:2:30;      %Array of SNR
snrx= 10.^(-snrdb/10);
ber4c1 = zeros(1,length(snrdb));
ber4c2 = zeros(1,length(snrdb));


bits=round(rand(1,l));  %Bits generation
bit = 0;                %Resultant bits

zqpsk = exp((1i*(pi/2)*[0:3]));          %qpsk "Unit Power"

M1 = length(zqpsk);

nbits1 = log2(M1);
sym1=zeros(1,l/nbits1);   %Symbol sent -> QPSK

for len = nbits1:nbits1:l
    switch (bi2de([bits(len) bits(len-1)]))
        case 0
            sym1(len/nbits1) = zqpsk(1);
        case 1
            sym1(len/nbits1) = zqpsk(2);         
        case 2
            sym1(len/nbits1) = zqpsk(3);            
        case 3
             sym1(len/nbits1) = zqpsk(4);           
    end
end

out = zeros(1,l/nbits1);     %Received symbols array

SD = 1/sqrt(2);
loop_Cond = sym1;

phase = exp(1i*(0));
symlen = size(sym1,2);

for Case=2
    for n=1:length(snrdb)
        H1l=(normrnd(0,SD,[1 symlen])+1i*normrnd(0,SD,[1 symlen]));
        H2l=(normrnd(0,SD,[1 symlen])+1i*normrnd(0,SD,[1 symlen]));
        H3l=(normrnd(0,SD,[1 symlen])+1i*normrnd(0,SD,[1 symlen]));
        H4l=(normrnd(0,SD,[1 symlen])+1i*normrnd(0,SD,[1 symlen]));            

        w1_00l = (sqrt(snrx(n))/2)*(randn(1,symlen) + 1i*randn(1,symlen));        %Noise Profile RH 0
        w1_10l = (sqrt(snrx(n))/2)*(randn(1,symlen) + 1i*randn(1,symlen));        %Noise Profile RV 0
        w1_01l = (sqrt(snrx(n))/2)*(randn(1,symlen) + 1i*randn(1,symlen));        %Noise Profile RH 1
        w1_11l = (sqrt(snrx(n))/2)*(randn(1,symlen) + 1i*randn(1,symlen));        %Noise Profile RV 1
        
        for loop=1:4:symlen        %Two symbols at every iteration.

            H1=H1l(loop);
            H2=H2l(loop);
            H3=H3l(loop);
            H4=H4l(loop);            
            
            w1_00 = w1_00l(loop);        %Noise Profile RH 0
            w1_10 = w1_10l(loop);        %Noise Profile RV 0
            w1_01 = w1_01l(loop);        %Noise Profile RH 1
            w1_11 = w1_11l(loop);        %Noise Profile RV 1

            q1=sym1(loop);                         % First symbol
            q2=sym1(loop+1);                 % Second symbol
            q3=sym1(loop+2);                 % Third symbol
            q4=sym1(loop+3);                       % fourth symbol

            q1o=sym1(loop)*phase;                         % First symbol rotated
            q2o=sym1(loop+1)*phase;                       % First symbol rotated
            q3o=sym1(loop+2)*phase;                       % First symbol rotated
            q4o=sym1(loop+3)*phase;                       % First symbol rotated
            
            if (Case == 1)
                RH0 = q1*H1 + q4*H3 - q2*conj(H2) - q3*conj(H4) + w1_00;      %Received at RH0            
                RV0 = q1*H2 + q4*H4 + q2*conj(H1) + q3*conj(H3) + w1_10;      %Received at RV0
                RH1 = conj(q2)*H1 - conj(q3)*H3 + conj(q1)*conj(H2) - conj(q4)*conj(H4) + w1_01;      %Received at RH1            
                RV1 = conj(q2)*H2 - conj(q3)*H4 - conj(q1)*conj(H1) + conj(q4)*conj(H3) + w1_11;      %Received at RV1
            elseif (Case == 2)
                RH0 = q1*H1 + q4*H3 - (q2o*conj(H2) + q3o*conj(H4)) + w1_00;      %Received at RH0            
                RV0 = q1*H2 + q4*H4 + (q2o*conj(H1) + q3o*conj(H3)) + w1_10;      %Received at RV0
                RH1 = conj(q2)*H1 - conj(q3)*H3 + (conj(q1o)*conj(H2) - conj(q4o)*conj(H4)) + w1_01;      %Received at RH1            
                RV1 = conj(q2)*H2 - conj(q3)*H4 + (-conj(q1o)*conj(H1) + conj(q4o)*conj(H3)) + w1_11;      %Received at RV1

            end

    %%%%%%%% Received Signal at Decoder %%%%%%%%%%%

            lenz = size(zqpsk);
            const1=100000;
            const2=100000;
            const3=100000;
            const4=100000;
            
            for lx=1:lenz(2)
                C = zqpsk(lx);
                if Case == 1
                    DR1 = -1*real(conj(RH0)*C*H1 + conj(RV0)*C*H2 + conj(RH1)*conj(C)*conj(H2) - conj(RV1)*conj(C)*conj(H1));       %Decoding first Symbol
                    DR2 = -1*real(-conj(RH0)*C*conj(H2) + conj(RV0)*C*conj(H1) + conj(RH1)*conj(C)*H1 + conj(RV1)*conj(C)*H2);       %Decoding second Symbol
                    DR3 = -1*real(-conj(RH0)*C*conj(H4) + conj(RV0)*C*conj(H3) - conj(RH1)*conj(C)*H3 - conj(RV1)*conj(C)*H4);       %Decoding third Symbol
                    DR4 = -1*real(conj(RH0)*C*H3 + conj(RV0)*C*H4 - conj(RH1)*conj(C)*conj(H4) + conj(RV1)*conj(C)*conj(H3));       %Decoding fourth Symbol
                elseif Case == 2
                    DR1 = -1*real(conj(RH0)*C*H1 + conj(RV0)*C*H2 + conj(C*phase)*(conj(RH1)*conj(H2) - conj(RV1)*conj(H1)));       %Decoding first Symbol
                    DR2 = -1*real(C*phase*(-conj(RH0)*conj(H2) + conj(RV0)*conj(H1)) + conj(RH1)*conj(C)*H1 + conj(RV1)*conj(C)*H2);       %Decoding second Symbol
                    DR3 = -1*real(C*phase*(-conj(RH0)*conj(H4) + conj(RV0)*conj(H3)) - conj(RH1)*conj(C)*H3 - conj(RV1)*conj(C)*H4);       %Decoding third Symbol
                    DR4 = -1*real(conj(RH0)*C*H3 + conj(RV0)*C*H4 + conj(C*phase)*(-conj(RH1)*conj(H4) + conj(RV1)*conj(H3)));       %Decoding fourth Symbol
                end
                
                if (DR1 < const1)
                   out(loop) = C;
                   const1 = DR1;
                end

                if (DR2 < const2)
                   out(loop+1) = C;
                   const2 = DR2;
                end

                if (DR3 < const3)
                   out(loop+2) = C;
                   const3 = DR3;
                end
                
                if (DR4 < const4)
                   out(loop+3) = C;
                   const4 = DR4;
                end
            end
        end

        for k=1:1:l/nbits1
            switch out(k)
                case zqpsk(1)
                   bit(nbits1*k-(nbits1-1):nbits1*k)=[0 0];
                case zqpsk(2)
                   bit(nbits1*k-(nbits1-1):nbits1*k)=[0 1];                
                case zqpsk(3)
                   bit(nbits1*k-(nbits1-1):nbits1*k)=[1 0];                
                case zqpsk(4)
                   bit(nbits1*k-(nbits1-1):nbits1*k)=[1 1];
            end
        end
        
        if Case==1
            ber4c1(n)=(l-sum(bits==bit))/l;
        elseif Case==2
            ber4c2(n)=(l-sum(bits==bit))/l;
        end 
    end
end
    
semilogy(snrdb,ber4c1,'m-*','linewidth',2)
hold on
semilogy(snrdb,ber4c2,'r-*','linewidth',2)
grid on
legend('QPSK,QPSK','QPSK,QPSKr')
xlim([-4 30])
ylim([10^-8 10^0])
title('2x1 system BER vs SnR Curve');
xlabel(' SNR (dB)') % x-axis label
ylabel(' BER ') % y-axis label
toc 
