bw=BW*cos(white/LIMIT*pi/2.0); set boxwidth bw; white=white+wd
red = sprintf("#%02X%02X%02X", 128+white/2, white, white)
green = sprintf("#%02X%02X%02X", white, 128+white/2, white)
blue = sprintf("#%02X%02X%02X", white, white, 128+white/2)
rep
if(white<LIMIT) reread

