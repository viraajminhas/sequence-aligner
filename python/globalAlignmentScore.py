def global_score_table(str1: str, str2: str, match: float, mismatch: float, gap: float) -> list[list[float]]:
    """
    Returns a 2-D array holding dynamic programming scores for global alignment.

    Parameters:
    - str1: str
    - str2: str
    - match: float
    - mismatch: float
    - gap: float

    Returns:
    - list[list[float]]: scoring table
    """
    if len(str1) == 0 or len(str2) == 0:
        raise ValueError("Error: empty string given to global_score_table.")

    num_rows = len(str1) + 1
    num_cols = len(str2) + 1

    score_table: list[list[float]] = []
    for i in range(num_rows):
        row: list[float] = []
        for j in range(num_cols):
            row.append(0.0)
        score_table.append(row)

    for j in range(1, num_cols):
        score_table[0][j] = float(j) * (-gap)
    for i in range(1, num_rows):
        score_table[i][0] = float(i) * (-gap)

    for i in range(1, num_rows):
        for j in range(1, num_cols):
            up_value = score_table[i-1][j] - gap
            left_value = score_table[i][j-1] - gap
            if str1[i-1] == str2[j-1]:
                diagonal_weight = match
            else:
                diagonal_weight = -mismatch
            diag_value = score_table[i-1][j-1] + diagonal_weight
            score_table[i][j] = max(up_value, left_value, diag_value)

    return score_table


def main() -> None:
    pass


if __name__ == "__main__":
    main()
