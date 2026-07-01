def local_score_table(str1: str, str2: str, match: float, mismatch: float, gap: float) -> list[list[float]]:
    """
    Returns the local alignment dynamic programming score table for the two strings (entries are floored at 0).

    Parameters:
    - str1: str
    - str2: str
    - match: float
    - mismatch: float
    - gap: float

    Returns:
    - list[list[float]]: the local alignment score table
    """
    if len(str1) == 0 or len(str2) == 0:
        raise ValueError("Error: empty string given to local_score_table.")

    num_rows = len(str1) + 1
    num_cols = len(str2) + 1

    score_table: list[list[float]] = []
    for i in range(num_rows):
        row: list[float] = []
        for j in range(num_cols):
            row.append(0.0)
        score_table.append(row)

    # The first row and first column stay 0: a local alignment is free to
    # start anywhere, so there is no penalty for the prefixes we skip.
    for i in range(1, num_rows):
        for j in range(1, num_cols):
            up_value = score_table[i-1][j] - gap
            left_value = score_table[i][j-1] - gap
            if str1[i-1] == str2[j-1]:
                diagonal_weight = match
            else:
                diagonal_weight = -mismatch
            diag_value = score_table[i-1][j-1] + diagonal_weight
            # Floor at 0: this is the one line that turns Needleman-Wunsch
            # into Smith-Waterman. Whenever the best running score would go
            # negative, we reset to 0 and let a new local alignment begin here.
            score_table[i][j] = max(0.0, up_value, left_value, diag_value)

    return score_table


def main() -> None:
    pass


if __name__ == "__main__":
    main()
