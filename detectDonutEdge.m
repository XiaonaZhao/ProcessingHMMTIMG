function [insideRing, outsideRing] = detectDonutEdge(DonutROImask)

bitwiseMask = ~DonutROImask;
centralMask = findMinDomain(bitwiseMask);

insideRing = edge(centralMask ,'canny'); % 'canny' would enroll some others.
outsideRing = edge((centralMask | (~bitwiseMask)), 'canny');
end
