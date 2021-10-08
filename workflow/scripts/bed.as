table gene_conversion
"Windows with evidence of gene conversion"
(
string  chrom;		"Chromosome for original alignment"
uint    chromStart;	"Start position of original alignment"
uint    chromEnd;	"End position of original alignment"
string  name;		"Name with mismatch delta"
uint    score;		"NA"
char[1]  strand;		"strand NA"
uint    thickStart;	"End"
uint    thickEnd;	"Start"
uint  reserved;		"RGB color, blue = acceptor, orange = donor"
string  donorChrom;		"Chromosome for donor alignment"
uint    donorStart;	"Start position of donor alignment"
uint    donorEnd;	"End position of donor alignment"
uint    mismatches;		"mismatches at original alignment"
uint    donorMismatches;		"mismatches at donor alignment"
float    perID_by_all;		"Percent identity of alignment at original alignment location"
float    donor_perID_by_all;		"Percent identity of alignment at donor location"
)