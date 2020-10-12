function potential = potentialLine(scanRate, sampleRate, highPotential, lowPotential)
scanRate = scanRate/1000;
scanRate = -scanRate/sampleRate;
potential1 = (highPotential : scanRate : lowPotential)';
end1 = potential1(end);
potential2 = ((lowPotential*2 - end1 - scanRate) : (-scanRate) : highPotential)';
% potential1 = (-0.6 : r : -1.1)';
% end1 = potential1(end);
% potential2 = ((-2.2 - end1 - r) : (-r) : -0.6)';

potential = [potential1' potential2' potential1' potential2']'; % 2c
% potential = [potential1' potential2']';% 1c
end