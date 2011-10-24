"""Get per-contig read counts from 454ReadStatus.txt produced by Newbler assembler.
We count 5' and 3' ends as 1/2 each to handle split reads."""

from MGT.Asm import *

contigReadCount454(inp=sys.stdin,out=sys.stdout)


