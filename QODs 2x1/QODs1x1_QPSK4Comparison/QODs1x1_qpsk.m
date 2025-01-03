clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The QPSK modulation %%%  

tic;        %Timer Started

%Two modulation Schemes QPSK & QPSKr
% QPSK <--- Z1
% QPSKr<--- Z2

l= floor(1000000/8)*8;            %No. of Bits @Note: should be divisble by nbits

snrdb=-4:2:30;      %Array of SNR
snrx= 10.^(-snrdb/20);
ber4c1 = zeros(1,length(snrdb));
ber4c2 = zeros(1,length(snrdb));


bits=round(rand(1,l));  %Bits generation
bit = 0;                %Resultant bits

zqpsk = exp((1i*(pi/2)*[0:3])+(0));          %qpsk "Unit Power"

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

for Case=1:2
    for n=1:length(snrdb)
        for loop=1:2:size(sym1,2)        %Two symbols at every iteration.

            Hx11=normrnd(0,SD,[1 1]);
            Hy11=normrnd(0,SD,[1 1]);
            H1=(Hx11+1i*Hy11);

            Hx22=normrnd(0,SD,[1 1]);
            Hy22=normrnd(0,SD,[1 1]);
            H2=(Hx22+1i*Hy22);

 
            w1_00 = snrx(n)*(1/sqrt(2))*(randn(1,1) + 1i*randn(1,1));        %Noise Profile RH 0
            w1_10 = snrx(n)*(1/sqrt(2))*(randn(1,1) + 1i*randn(1,1));        %Noise Profile RV 0

            q1=sym1(loop);                   % First symbol
            q2=sym1(loop+1);                 % Second symbol
            
            if (Case == 1)
                RH0 = q1*H1 - q2*conj(H2) + w1_00;      %Received at RH0            
                RV0 = q1*H2 + q2*conj(H1) + w1_10;      %Received at RV0
            elseif (Case == 2)
                RH0 = q1*H1 - q2*exp(1i*(pi/4))*conj(H2) + w1_00;      %Received at RH0            
                RV0 = q1*H2 + q2*exp(1i*(pi/4))*conj(H1) + w1_10;      %Received at RV0
            end

    %%%%%%%% Received Signal at Decoder %%%%%%%%%%%

            lenz = size(zqpsk);
            const1=100000;
            const2=100000;
            
            for lx=1:lenz(2)
                C = zqpsk(lx);
                
                if Case == 1
                    DR1 = -1*real(C*H1*conj(RH0) + conj(C)*RV0*conj(H2));       %Decoding first Symbol
                    DR2 = -1*real(-C*conj(RH0)*conj(H2) + conj(C)*RV0*H1);       %Decoding second Symbol
                elseif Case == 2
                    DR1 = -1*real(C*H1*conj(RH0) + conj(C)*RV0*conj(H2));                                      %Decoding first Symbol
                    DR2 = -1*real(-C*exp(1i*(pi/4))*conj(RH0)*conj(H2) + conj(C*exp(1i*(pi/4)))*RV0*H1);       %Decoding second Symbol
                end
                
                if (DR1 < const1)
                   out(loop) = C;
                   const1 = DR1;
                end

                if (DR2 < const2)
                   out(loop+1) = C;
                   const2 = DR2;
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
title('1x1 system BER vs SnR Curve');
xlabel(' SNR (dB)') % x-axis label
ylabel(' BER ') % y-axis label
toc 
