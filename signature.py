import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze

Analyze.cosmic_fit(samples="signature_matrix.txt", 
                   output="signature",
                   input_type="matrix",
                   genome_build="GRCh37",
                   cosmic_version=3.3)