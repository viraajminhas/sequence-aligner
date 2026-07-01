from localAlignmentScore import local_score_table

type Alignment = list[str]


def local_alignment(str1: str, str2: str, match: float, mismatch: float, gap: float) -> tuple[Alignment, int, int, int, int]:
    """
    Returns a highest-scoring local alignment of the two strings, with the start and end indices of the aligned substring in each string.

    Parameters:
    - str1: str
    - str2: str
    - match: float
    - mismatch: float
    - gap: float

    Returns:
    - tuple[Alignment, int, int, int, int]: the alignment, then start1, end1, start2, end2
    """
    alignment = ["", ""]

    scoring_matrix = local_score_table(str1, str2, match, mismatch, gap)

    # Find the highest-scoring cell anywhere in the table: that is where the
    # best local alignment ends.
    best_score = 0.0
    r = 0
    c = 0
    for i in range(len(scoring_matrix)):
        for j in range(len(scoring_matrix[i])):
            if scoring_matrix[i][j] > best_score:
                best_score = scoring_matrix[i][j]
                r = i
                c = j

    end1 = r
    end2 = c

    # Backtrack from the best cell until we hit a 0, which marks the start of
    # the local alignment.
    while r > 0 and c > 0 and scoring_matrix[r][c] > 0:
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
            raise ValueError("Error: unexpected state in local alignment backtracking.")

    start1 = r
    start2 = c

    return alignment, start1, end1, start2, end2


def main() -> None:
    pass


if __name__ == "__main__":
    main()
