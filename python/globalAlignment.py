from globalAlignmentScore import global_score_table

type Alignment = list[str]


def global_alignment(str1: str, str2: str, match: float, mismatch: float, gap: float) -> Alignment:
    """
    Returns a maximum score global alignment of two strings using the given penalties.

    Parameters:
    - str1: str
    - str2: str
    - match: float
    - mismatch: float
    - gap: float

    Returns:
    - Alignment: list of two aligned strings
    """
    alignment = ["", ""]

    scoring_matrix = global_score_table(str1, str2, match, mismatch, gap)

    r = len(str1)
    c = len(str2)

    while r > 0 and c > 0:
        if scoring_matrix[r][c] == scoring_matrix[r-1][c] - gap:
            # UP: symbol from string 1 against gap
            alignment[0] = str1[r-1] + alignment[0]
            alignment[1] = "-" + alignment[1]
            r -= 1
        elif scoring_matrix[r][c] == scoring_matrix[r][c-1] - gap:
            # LEFT: gap against symbol of string 2
            alignment[0] = "-" + alignment[0]
            alignment[1] = str2[c-1] + alignment[1]
            c -= 1
        elif str1[r-1] != str2[c-1] and scoring_matrix[r][c] == scoring_matrix[r-1][c-1] - mismatch:
            # mismatch case
            alignment[0] = str1[r-1] + alignment[0]
            alignment[1] = str2[c-1] + alignment[1]
            r -= 1
            c -= 1
        elif str1[r-1] == str2[c-1] and scoring_matrix[r][c] == scoring_matrix[r-1][c-1] + match:
            # match case
            alignment[0] = str1[r-1] + alignment[0]
            alignment[1] = str2[c-1] + alignment[1]
            r -= 1
            c -= 1
        else:
            raise ValueError("Error: unexpected state in global alignment backtracking.")

    while r > 0:
        # move up until origin: align symbol from string 1 against gap
        alignment[0] = str1[r-1] + alignment[0]
        alignment[1] = "-" + alignment[1]
        r -= 1

    while c > 0:
        # move left until origin: gap against symbol of string 2
        alignment[0] = "-" + alignment[0]
        alignment[1] = str2[c-1] + alignment[1]
        c -= 1

    return alignment


def main() -> None:
    pass


if __name__ == "__main__":
    main()
