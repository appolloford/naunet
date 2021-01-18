def request_user_check(question: str, exit_option: str) -> None:
    answer = input(question)
    if answer.lower() == exit_option.lower():
        exit
