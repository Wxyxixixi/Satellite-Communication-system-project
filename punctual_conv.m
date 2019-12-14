% 
clear all;
puncturepattern12=[1;1];%puncture pattern when the code rate=1/2
puncturepattern23=[1;1;0;1];                                                                                                   %2/3
puncturepattern56=[1;1;0;1;1;0;0;1;1;0];  
puncturepatternSelect=[puncturepattern12;puncturepattern23;puncturepattern56];
PC=[1,2;3,6;7,16];
tracebackdepthSelect=[60,75,120];
coderateSelect=[1/2,2/3,5/6];
M=4;
EbN0=-4:0.1:16;
BER_QPSK=zeros(3,length(EbN0));
EsNo=EbN0+10*log10(log2(M));
for b=1:3
EbNoratesemilogy=zeros(1,length(EbN0));  
hConvEnc = comm.ConvolutionalEncoder(poly2trellis(7,[133 171]));       %maybe need to change the 133 171
hConvEnc.PuncturePatternSource = 'Property';
hConvEnc.PuncturePattern =puncturepatternSelect(PC(b,1):PC(b,2),1);
hMod = comm.PSKModulator(M, 'BitInput',true); 

hDemod = comm.PSKDemodulator(M, 'BitOutput',true);
hVitDec = comm.ViterbiDecoder(poly2trellis(7,[133 171]),...
    'InputFormat','Hard');

hVitDec.PuncturePatternSource  = 'Property';                           %property not sure
hVitDec.PuncturePattern = hConvEnc.PuncturePattern;
hVitDec.TracebackDepth = tracebackdepthSelect(b);

hErrorCalc = comm.ErrorRate('ReceiveDelay',hVitDec.TracebackDepth);


EbNoEncoderOutput = EbN0 +10*log10(coderateSelect(b));

frameLength =45000;                                                     %framelength/maxnumtransmissions
targetErrors = 100; %100 errors max
maxNumTransmissions =10000000;
BERVec = zeros(3,length(EbNoEncoderOutput));
%BERp=zeros(3,length(EbNoEncoderOutput));

for n=1:length(EbNoEncoderOutput)
    reset(hErrorCalc);
    reset(hConvEnc);
    reset(hVitDec);
    [EbNoEncoderOutput(n),coderateSelect(b)]
    hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Es/No)','EsNo',EsNo(n));           %signalpower not sure

    while(BERVec(2,n) < targetErrors)&&(BERVec(3,n)<maxNumTransmissions)
        data = randi([0 1],frameLength,1);                               %use framelength to create random number?
        encData = step(hConvEnc, data);
        modData = step(hMod,encData);
        channelOutput = step(hChan,modData);
        DemodData = step(hDemod, channelOutput);
        decData = step(hVitDec, (DemodData));
        %d = tracebackdepthSelect;
        BERVec(:,n) = step(hErrorCalc, data, decData);
   
    end
    if BERVec(1,n)<2e-6
        
       BER_QPSK(b,:) = BERVec(1,:);
%BERp(b,:)=BERVec(1,:);
save('BPSK_sim_3.mat','EbN0','BER_QPSK');
load('BPSK_sim_3.mat');

figure(1)
semilogy(EbN0,BER_QPSK(b,:),'-');
hold on;

grid on;
title('Part2 question1 BER verus Eb/N0 curve without RS outer code ');
ylabel('BER');
xlabel('EbN0 Ratio'); 
       
        
    break
    end
end

%BERp(b,:)=BERVec(1,:);
save('BPSK_sim_3.mat','EbN0','BER_QPSK');
load('BPSK_sim_3.mat');
 
end
puncturepattern12=[1;1];                                 %puncture pattern when the code rate=1/2
puncturepattern23=[1;1;0;1];                                                                                                   %2/3
puncturepattern56=[1;1;0;1;1;0;0;1;1;0];                                                                                %5/6
puncturepatternSelect=[puncturepattern12;puncturepattern23;puncturepattern56]; % a column array (2+4+10)*1 
PC=[1,2;3,6;7,16];                                                     
coderateSelect=[1/2,2/3,5/6];
M=4;
EbN0=-4:0.5:16;
BER8=zeros(3,length(EbN0));
for b=1:3
  EbNoratesemilogy=zeros(1,length(EbN0)); 
 hErrorCalc = comm.ErrorRate(); 
%hError = comm.ErrorRate();

 hMod = comm.PSKModulator(M, 'BitInput',true); 
 hDemod = comm.PSKDemodulator(M, 'BitOutput',true);
%hQPSKMod=comm.QPSKModulator;
%hQPSKDemod=comm.QPSKDemodulator;

hChan= comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)');
hChan.BitsPerSymbol=2;

BERVec1=zeros(3,length(EbN0));
EbNoEncoderOutput = EbN0 +10*log10(coderateSelect(b));
for n=1:length(EbNoEncoderOutput)
     reset(hErrorCalc);
     [EbNoEncoderOutput(n),coderateSelect(b)]
    hChan.EbNo = EbNoEncoderOutput(n);
    
     %hChanQPSK.EbNo =SNR(n);
     
     while (BERVec1(2,n)<100)&&(BERVec1(3,n)<5000000)
         data0 = randi([0,1],64800, 1);
         
         modQ=step(hMod,data0);
         chanQ=step(hChan,modQ);
         demodQ=step(hDemod,chanQ);
        

         BERVec1(:,n)= step(hErrorCalc, data0, demodQ);
         BER8(b,n)=BERVec1(1,n);
     end
end
semilogy(EbN0, BER8(b,:),'-');


hold on;
end
 legend('1/2puncturingrate','2/3puncturingrate','5/6puncturingrate','1/2uncodedQPSK','2/3uncodedQPSK', '5/6uncodedQPSK');
 
 
 
 
 
 
 

 