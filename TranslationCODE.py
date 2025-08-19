
from typing import Literal, Tuple, Listvhfccccccccccccccccccccccccccccccccccccc


# Standard RNA codon table → one-letter amino acids
CODON_TO_AA = {
    # Phe / Leu
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    # Ile / Met
    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
    # Val
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    # Ser / Pro / Thr / Ala
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S","AGU":"S","AGC":"S",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    # Tyr / His / Gln / Asn / Lys / Asp / Glu / Cys / Trp
    "UAU":"Y","UAC":"Y",
    "CAU":"H","CAC":"H",
    "CAA":"Q","CAG":"Q",
    "AAU":"N","AAC":"N",
    "AAA":"K","AAG":"K",
    "GAU":"D","GAC":"D",
    "GAA":"E","GAG":"E",
    "UGU":"C","UGC":"C",
    "UGG":"W",
    # Arg / Gly
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G",
    # Stops
    "UAA":"*","UAG":"*","UGA":"*",
}

AA1_TO_AA3 = {
    "A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys","Q":"Gln","E":"Glu",
    "G":"Gly","H":"His","I":"Ile","L":"Leu","K":"Lys","M":"Met","F":"Phe",
    "P":"Pro","S":"Ser","T":"Thr","W":"Trp","Y":"Tyr","V":"Val","*":"Stop","X":"Xaa"
}

def translate_mrna(seq: str,
                   frame: int = 0,
                   *,
                   start_at_first_aug: bool = False,
                   to_stop: bool = False,
                   alphabet: Literal["1", "3"] = "1",
                   sep: str | None = None,
                   unknown: str = "X",
                   strict: bool = False) -> Tuple[str, int]:
    """
    Translate an mRNA (or DNA) sequence to protein using the standard genetic code.

    Parameters
    ----------
    seq : str
        Nucleotide sequence (RNA or DNA). Whitespace is ignored; case-insensitive.
    frame : int
        Reading frame (0, 1, or 2) counted from the first base of `seq`.
    start_at_first_aug : bool
        If True, scan codon-wise in the chosen frame and start translation at the first AUG.
    to_stop : bool
        If True, stop translation at the first stop codon (UAA/UAG/UGA) without including '*'.
    alphabet : {"1","3"}
        "1" for one-letter amino acids (default), "3" for three-letter (e.g., Met-Ser-...).
    sep : str | None
        Separator when alphabet == "3". Default: "-" for 3-letter, "" for 1-letter.
    unknown : str
        Symbol to use when a codon contains ambiguous base(s) (e.g., N, R, Y). Default "X".
    strict : bool
        If True, raise ValueError if trailing bases are not a full codon; otherwise they're ignored.

    Returns
    -------
    protein : str
        The translated protein.
    start_index : int
        Zero-based nucleotide index in the ORIGINAL input where translation began.

    Notes
    -----
    - DNA 'T' is converted to RNA 'U'.
    - Ambiguous codons (non-ACGU) yield `unknown`.
    """
    if frame not in (0,1,2):
        raise ValueError("frame must be 0, 1, or 2")

    # Normalize input
    rna = "".join(seq.upper().split()).replace("T", "U")
    n = len(rna)

    # Decide where to start
    start = frame
    if start_at_first_aug:
        found = False
        for i in range(frame, n - 2, 3):
            if rna[i:i+3] == "AUG":
                start = i
                found = True
                break
        if not found:
            # No AUG in frame: fall back to provided frame start
            start = frame

    # Handle trailing bases
    remainder = (n - start) % 3
    if remainder and strict:
        raise ValueError(f"Sequence length from start {start} not divisible by 3.")
    limit = n - remainder

    aa1_list: List[str] = []
    stop_hit = False

    for i in range(start, limit, 3):
        codon = rna[i:i+3]
        if any(b not in "ACGU" for b in codon):
            aa = unknown
        else:
            aa = CODON_TO_AA.get(codon, unknown)

        if aa == "*":
            if to_stop:
                stop_hit = True
                break
            else:
                aa1_list.append("*")
        else:
            aa1_list.append(aa)

    # Format output
    if alphabet == "1":
        protein = "".join(aa1_list)
    else:
        sep = "-" if sep is None else sep
        protein = sep.join(AA1_TO_AA3.get(a, AA1_TO_AA3["X"]) for a in aa1_list)

    return protein, start
