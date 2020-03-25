# Example PBS scripts

Parallel multiscaling is an extremely memory-intensive exercise. On Raijin, multiscaling jobs used to be run on the 
`megamem` queue (see example PBS script for Raijin), which had ~3TB of RAM per node (32 cores). 

Unfortunately, there is no such equivalent on Gadi and we use a simple workaround: we acquire 2 `hugemem` nodes 
(each with 48 cores and 1.5TB RAM), but run only 32 processes on these two nodes (16 per node, as directed by
`--map-by ppr:16:node` in the example pbs script) -- this gives us the same RAM per core as that on Raijin.
