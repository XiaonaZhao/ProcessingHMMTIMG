function Voltage  = calculateVolt(Current, BeginVolt, EndVolt)

Dots = length(Current);
Voltage1 = BeginVolt;
Voltage2 = EndVolt;
VoltDots = (linspace(Voltage1,Voltage2,...
    (Dots + mod(Dots, 2))/2 ))';% linspace is row vector
if mod(Dots, 2)
    Voltage = [VoltDots; flipud(VoltDots(1:end-1))];
else
    Voltage = [VoltDots; flipud(VoltDots)];
end
end