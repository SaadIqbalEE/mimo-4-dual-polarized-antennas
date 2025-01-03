clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The QPSK modulation %%%  

tic;        %Timer Started

l=100000;            %No. of Bits @Note: should be divisble by nbits
snrdb=-4:2:30;      %Array of SNR
snrx= 10.^(-snrdb/20);
ber = zeros(1,length(snrdb));
bits=round(rand(1,l));  %Bits generation
bit = 0;                %Resultant bits

zqpsk = (1/sqrt(2))*[1+1i -1+1i 1-1i -1-1i];    %Grey coding "Unit Power"
M = length(zqpsk);

nbits = log2(M);
sym=zeros(1,l/nbits);   %Symbol sent

for len = nbits:nbits:l
    switch (bi2de([bits(len) bits(len-1)]))
        case 0
            sym(len/nbits) = zqpsk(1);
        case 1
            sym(len/nbits) = zqpsk(2);         
        case 2
            sym(len/nbits) = zqpsk(3);            
        case 3
             sym(len/nbits) = zqpsk(4);           
    end
end

out = zeros(1,l/nbits);     %Received symbols array
SD = 1/sqrt(2);

for n=1:length(snrdb)
    for loop=1:2:size(sym,2)        %Two symbols at every iteration.
       
        Hx11=normrnd(0,SD,[1 1]);
        Hy11=normrnd(0,SD,[1 1]);
        H1=(Hx11+1i*Hy11);

        Hx22=normrnd(0,SD,[1 1]);
        Hy22=normrnd(0,SD,[1 1]);
        H2=(Hx22+1i*Hy22);
        
        w1_0 = snrx(n)*(1/sqrt(2))*(randn(1,1) + 1i*randn(1,1));        %Noise Profile RH
        w1_1 = snrx(n)*(1/sqrt(2))*(randn(1,1) + 1i*randn(1,1));        %Noise Profile RV
        
        q1=sym(loop);                   % First symbol
        q2=sym(loop+1);                 % Second symbol

        RH = q1*H1 - q2*conj(H2) + w1_0;      %Received at RH            
        RV = q1*H2 + q2*conj(H1) + w1_1;      %Received at RV

%%%%%%%% Received Signal at Decoder %%%%%%%%%%%
        
        const=-100000;
        for z1=1:4
            for z2=1:4
                
                C1 = zqpsk(z1);
                C2 = zqpsk(z2);
                
                DR = (real(conj(RH)*(C1*H1-C2*conj(H2)) + RV*(conj(C1)*conj(H2)+conj(C2)*H1)));       %Decoding Symbol
                if (DR > const)
                    temp = [z1 z2];
                    const= DR;
                end
            end
        end
        out(loop) = zqpsk(temp(1));
        out(loop+1) = zqpsk(temp(2));
        
    end

    for k=1:1:l/nbits
        switch out(k)
            case zqpsk(1)
               bit(nbits*k-(nbits-1):nbits*k)=[0 0];
            case zqpsk(2)
               bit(nbits*k-(nbits-1):nbits*k)=[0 1];                
            case zqpsk(3)
               bit(nbits*k-(nbits-1):nbits*k)=[1 0];                
            case zqpsk(4)
               bit(nbits*k-(nbits-1):nbits*k)=[1 1];
        end
    end
    ber(n)=(l-sum(bits==bit))/l;
end

semilogy(snrdb,ber,'m-*','linewidth',2)
grid on
xlim([-4 30])
ylim([10^-8 10^0])
title('1x1 system QPSK BER vs SnR Curve');
xlabel(' SNR (dB)') % x-axis label
ylabel(' BER ') % y-axis label
toc 


        