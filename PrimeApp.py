import os
import tkinter as tk
from datetime import datetime
from tkinter import messagebox, simpledialog
import time
import matplotlib.pyplot as plt

class PrimeCalculator:
    def __init__(self):
        self.previous_results = []

    def sieve_of_eratosthenes(self, start, end):
        primes = []
        sieve = [True] * (end + 1)
        sieve[0] = False
        if end > 0:
            sieve[1] = False
        for p in range(2, end + 1):
            if sieve[p]:
                if p >= start:
                    primes.append(p)
                for i in range(p * p, end + 1, p):
                    sieve[i] = False
        return primes

    def sieve_of_atkin(self, start, end):
        sieve = [False] * (end + 1)
        primes = []
        if end > 2:
            primes.append(2)
        if end > 3:
            primes.append(3)
        for x in range(1, int(end ** 0.5) + 1):
            for y in range(1, int(end ** 0.5) + 1):
                n = 4 * x**2 + y**2
                if n <= end and (n % 12 == 1 or n % 12 == 5):
                    sieve[n] = not sieve[n]
                n = 3 * x**2 + y**2
                if n <= end and n % 12 == 7:
                    sieve[n] = not sieve[n]
                n = 3 * x**2 - y**2
                if x > y and n <= end and n % 12 == 11:
                    sieve[n] = not sieve[n]
        for n in range(5, int(end ** 0.5) + 1):
            if sieve[n]:
                for k in range(n**2, end + 1, n**2):
                    sieve[k] = False
        for n in range(5, end + 1):
            if sieve[n] and n >= start:
                primes.append(n)
        return primes

    def is_prime(self, n):
        if n <= 1:
            return False
        if n <= 3:
            return True
        if n % 2 == 0 or n % 3 == 0:
            return False
        i = 5
        while i * i <= n:
            if n % i == 0 or n % (i + 2) == 0:
                return False
            i += 6
        return True

    def find_primes_simple(self, start, end):
        primes = []
        for num in range(start, end + 1):
            if self.is_prime(num):
                primes.append(num)
        return primes

    def segmented_sieve(self, l, r):
        import math

        def simple_sieve(limit):
            mark = [True] * (limit + 1)
            primes = []
            for p in range(2, limit + 1):
                if mark[p]:
                    primes.append(p)
                    for i in range(p * p, limit + 1, p):
                        mark[i] = False
            return primes

        limit = int(math.sqrt(r)) + 1
        primes = simple_sieve(limit)

        mark = [True] * (r - l + 1)

        for p in primes:
            low_limit = max(p * p, (l + p - 1) // p * p)
            for j in range(low_limit, r + 1, p):
                mark[j - l] = False

        prime_numbers = [num for num, is_prime in zip(range(l, r + 1), mark) if is_prime and num > 1]

        return prime_numbers

    def trial_division(self, start, end):
        primes = []
        for n in range(start, end + 1):
            if n <= 1:
                continue
            prime = True
            for d in range(2, int(n ** 0.5) + 1):
                if n % d == 0:
                    prime = False
                    break
            if prime:
                primes.append(n)
        return primes

    def wheel_factorization(self, start, end):
        def gcd(a, b):
            while b:
                a, b = b, a % b
            return a

        def lcm(a, b):
            return abs(a * b) // gcd(a, b)

        wheel = [2, 3, 5]
        w_lcm = 1
        for x in wheel:
            w_lcm = lcm(w_lcm, x)

        steps = [x for x in range(w_lcm) if all(x % p for p in wheel)]
        primes = [x for x in wheel if x >= start]
        is_prime = [True] * (end + 1)

        for x in range(2, int(end ** 0.5) + 1):
            if is_prime[x]:
                for i in range(x * x, end + 1, x):
                    is_prime[i] = False

        for x in range(max(start, wheel[-1] + 1), end + 1):
            if is_prime[x]:
                primes.append(x)

        return primes

    def count_primes(self, start, end, methods_selected):
        all_primes = set()
        results = []
        for method in methods_selected:
            start_time = time.perf_counter()
            if method == "Решето Эратосфена":
                primes = self.sieve_of_eratosthenes(start, end)
            elif method == "Решето Аткина":
                primes = self.sieve_of_atkin(start, end)
            elif method == "Простой перебор":
                primes = self.find_primes_simple(start, end)
            elif method == "Сегментированное решето":
                primes = self.segmented_sieve(start, end)
            elif method == "Пробное деление":
                primes = self.trial_division(start, end)
            elif method == "Факторизация колеса":
                primes = self.wheel_factorization(start, end)
            else:
                continue
            elapsed_time = time.perf_counter() - start_time
            results.append((method, len(primes), elapsed_time))
            self.previous_results.append((method, len(primes), elapsed_time))
            all_primes.update(primes)
        return list(all_primes), results

    def find_primes_simple(self, start, end, step, method):
        all_primes = []
        for i in range(start, end + 1, step):
            primes = []
            start_time = time.perf_counter()
            if method == "Решето Эратосфена":
                primes = self.sieve_of_eratosthenes(i, min(i + step - 1, end))
            elif method == "Решето Аткина":
                primes = self.sieve_of_atkin(i, min(i + step - 1, end))
            elif method == "Простой перебор":
                primes = self.find_primes_simple(i, min(i + step - 1, end))
            elif method == "Сегментированное решето":
                primes = self.segmented_sieve(i, min(i + step - 1, end))
            elif method == "Пробное деление":
                primes = self.trial_division(i, min(i + step - 1, end))
            elif method == "Факторизация колеса":
                primes = self.wheel_factorization(i, min(i + step - 1, end))
            elapsed_time = time.perf_counter() - start_time
            self.previous_results.append((method, len(primes), elapsed_time))
            all_primes.extend(primes)
        return all_primes

    def write_to_file(self, primes, filename):
        with open(filename, "w") as file:
            for method, count, elapsed_time in self.previous_results:
                file.write(f"{elapsed_time} {count}\n")

    def write_to_file1(self, primes,  filename):

        with open(filename, "w") as file:
            file.write("\n".join(map(str, sorted(primes))))

    def analyze_results(self):
        methods = list(set(method for method, _, _ in self.previous_results))
        data = {method: {"counts": [], "times": []} for method in methods}

        for method, count, elapsed_time in self.previous_results:
            data[method]["counts"].append(count)
            data[method]["times"].append(elapsed_time)

        plt.figure(figsize=(10, 6))
        for method in methods:
            plt.plot(data[method]["counts"], data[method]["times"], label=method, marker='o')

        plt.xlabel('Количество простых чисел')
        plt.ylabel('Время выполнения (секунды)')
        plt.title('Анализ времени выполнения различных алгоритмов')
        plt.legend()
        plt.grid(True)
        plt.show()

class PrimeApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Выбор программы")
        self.root.geometry("300x200")
        self.root.resizable(width=False, height=False)

        self.label_choose = tk.Label(self.root, text="Выберите программу:")
        self.label_choose.pack(pady=10)

        self.button_normal = tk.Button(self.root, text="Обычный поиск", command=self.open_normal_search)
        self.button_normal.pack(pady=5)

        self.button_step = tk.Button(self.root, text="Поиск с шагом", command=self.open_step_search)
        self.button_step.pack(pady=5)

        self.button_exit = tk.Button(self.root, text="Выход", command=self.root.destroy)
        self.button_exit.pack(pady=5)

    def open_normal_search(self):
        self.root.destroy()
        root_normal = tk.Tk()
        app_normal = NormalPrimeApp(root_normal)
        root_normal.mainloop()

    def open_step_search(self):
        self.root.destroy()
        root_step = tk.Tk()
        app_step = StepPrimeApp(root_step)
        root_step.mainloop()

class NormalPrimeApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Обычный поиск простых чисел")
        self.root.geometry("500x700")
        self.root.resizable(width=False, height=False)

        self.prime_calculator = PrimeCalculator()

        self.label_start = tk.Label(self.root, text="Начальное число:")
        self.label_start.pack(pady=10)

        self.entry_start = tk.Entry(self.root)
        self.entry_start.pack(pady=10)

        self.label_end = tk.Label(self.root, text="Конечное число:")
        self.label_end.pack(pady=10)

        self.entry_end = tk.Entry(self.root)
        self.entry_end.pack(pady=10)

        self.label_method = tk.Label(self.root, text="Выберите методы:")
        self.label_method.pack(pady=10)

        self.methods = [
            "Решето Эратосфена",
            "Решето Аткина",
            "Простой перебор",
            "Сегментированное решето",
            "Пробное деление",
            "Факторизация колеса"
        ]
        self.check_vars = [tk.BooleanVar() for _ in self.methods]

        for var, method in zip(self.check_vars, self.methods):
            frame = tk.Frame(self.root)
            check = tk.Checkbutton(frame, text=method, variable=var)
            check.pack(side=tk.LEFT)
            info_button = tk.Button(frame, text="?", command=lambda m=method: self.show_method_info(m))
            info_button.pack(side=tk.LEFT, padx=5)
            frame.pack(pady=5)

        self.button_recommendation = tk.Button(self.root, text="Рекомендации", command=self.show_algorithm_recommendation)
        self.button_recommendation.pack(pady=10)

        self.button_description = tk.Button(self.root, text="Описание всех алгоритмов", command=self.show_algorithm_description)
        self.button_description.pack(pady=10)

        self.button_count = tk.Button(self.root, text="Подсчет", command=self.count_primes)
        self.button_count.pack(pady=10)

        self.button_analyze = tk.Button(self.root, text="Анализ результатов", command=self.prime_calculator.analyze_results)
        self.button_analyze.pack(pady=10)

        self.button_exit = tk.Button(self.root, text="Выход", command=self.root.destroy)
        self.button_exit.pack(pady=10)

    def show_algorithm_recommendation(self):
        recommendation = (
            "Для небольших чисел рекомендуется использовать методы простого перебора и пробного деления. "
            "Для больших чисел лучше всего подходят решето Эратосфена и сегментированное решето."
        )
        messagebox.showinfo("Рекомендации по алгоритмам", recommendation)

    def show_method_info(self, method):
        info = {
            "Решето Эратосфена": (
                "Решето Эратосфена:\n"
                "Классический алгоритм поиска простых чисел. Эффективен для больших чисел.\n"
                "Временная сложность: O(n log log n)\n"
                "Рекомендация: Используйте для больших чисел."
            ),
            "Решето Аткина": (
                "Решето Аткина:\n"
                "Оптимизированный алгоритм поиска простых чисел. Эффективен для больших чисел.\n"
                "Временная сложность: O(n / log log n)\n"
                "Рекомендация: Используйте для больших чисел."
            ),
            "Простой перебор": (
                "Простой перебор:\n"
                "Наивный метод проверки каждого числа на простоту. Эффективен для небольших чисел.\n"
                "Временная сложность: O(n^2)\n"
                "Рекомендация: Используйте для небольших чисел."
            ),
            "Сегментированное решето": (
                "Сегментированное решето:\n"
                "Оптимизированный метод, использующий решето Эратосфена для сегментов диапазона.\n"
                "Временная сложность: O(n log log n)\n"
                "Рекомендация: Используйте для очень больших чисел."
            ),
            "Пробное деление": (
                "Пробное деление:\n"
                "Метод проверки делимости числа на все меньшие простые числа.\n"
                "Временная сложность: O(n^1.5)\n"
                "Рекомендация: Используйте для небольших чисел."
            ),
            "Факторизация колеса": (
                "Факторизация колеса:\n"
                "Оптимизированный метод пробного деления, исключающий делители, которые являются известными простыми числами.\n"
                "Временная сложность: O(n^1.5)\n"
                "Рекомендация: Используйте для чисел среднего размера."
            ),
        }
        messagebox.showinfo(method, info.get(method, "Информация не найдена."))

    def show_algorithm_description(self):
        descriptions = (
            "Решето Эратосфена:\n"
            "Классический алгоритм поиска простых чисел. Эффективен для больших чисел.\n"
            "Временная сложность: O(n log log n)\n\n"
            "Решето Аткина:\n"
            "Оптимизированный алгоритм поиска простых чисел. Эффективен для больших чисел.\n"
            "Временная сложность: O(n / log log n)\n\n"
            "Простой перебор:\n"
            "Наивный метод проверки каждого числа на простоту. Эффективен для небольших чисел.\n"
            "Временная сложность: O(n^2)\n\n"
            "Сегментированное решето:\n"
            "Оптимизированный метод, использующий решето Эратосфена для сегментов диапазона.\n"
            "Временная сложность: O(n log log n)\n\n"
            "Пробное деление:\n"
            "Метод проверки делимости числа на все меньшие простые числа.\n"
            "Временная сложность: O(n^1.5)\n\n"
            "Факторизация колеса:\n"
            "Оптимизированный метод пробного деления, исключающий делители, которые являются известными простыми числами.\n"
            "Временная сложность: O(n^1.5)\n"
        )
        messagebox.showinfo("Описание алгоритмов", descriptions)

    def count_primes(self):
        try:
            start = int(self.entry_start.get())
            end = int(self.entry_end.get())
        except ValueError:
            messagebox.showerror("Ошибка", "Пожалуйста, введите корректные числа")
            return

        methods_selected = [method for method, var in zip(self.methods, self.check_vars) if var.get()]

        if not methods_selected:
            messagebox.showerror("Ошибка", "Пожалуйста, выберите хотя бы один метод")
            return

        primes, results = self.prime_calculator.count_primes(start, end, methods_selected)
        self.prime_calculator.write_to_file1(primes, "output.txt")

        result_str = "\n\n".join(f"Метод: {method}\nКоличество простых чисел: {count}\nВремя выполнения: {time:.6f} сек"
                                 for method, count, time in results)
        messagebox.showinfo("Результаты", result_str)

class AnalysisApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Анализ результатов")
        self.root.geometry("500x500")
        self.root.resizable(width=False, height=False)

        self.selected_files = []

        self.label_choose = tk.Label(self.root, text="Выберите файлы для анализа:")
        self.label_choose.pack(pady=10)

        self.file_button = tk.Button(self.root, text="Выбрать файл", command=self.select_file)
        self.file_button.pack(pady=10)

        self.analyze_button = tk.Button(self.root, text="Анализировать", command=self.analyze_files)
        self.analyze_button.pack(pady=10)

        self.selected_files_label = tk.Label(self.root, text="")
        self.selected_files_label.pack(pady=10)

        self.result_label = tk.Label(self.root, text="")
        self.result_label.pack(pady=10)

        self.back_button = tk.Button(self.root, text="Назад", command=self.reset)
        self.back_button.pack(pady=10)

        self.button_exit = tk.Button(self.root, text="Выход", command=self.root.destroy)
        self.button_exit.pack(pady=10)

        self.figure = None  # To store the plot figure

    def select_file(self):
        try:
            filename = simpledialog.askstring("Input", "Введите имя файла для анализа:", parent=self.root)
            if filename:
                filename = filename.strip()
                if os.path.isfile(filename):
                    self.selected_files.append(filename)
                    self.selected_files_label.config(text="Выбранные файлы:\n" + "\n".join([os.path.basename(f) for f in self.selected_files]))
                else:
                    messagebox.showerror("Ошибка", f"Файл не найден: {filename}")
            else:
                self.selected_files_label.config(text="")
        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка при выборе файла: {e}")

    def analyze_files(self):
        if not self.selected_files:
            messagebox.showerror("Ошибка", "Пожалуйста, выберите файлы для анализа")
            return

        try:
            self.figure, ax = plt.subplots(figsize=(10, 6))
            colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
            for i, selected_file in enumerate(self.selected_files):
                with open(selected_file, "r") as file:
                    data = file.readlines()
                    counts = []
                    times = []
                    accumulated_count = 0
                    accumulated_time = 0
                    for line in reversed(data):
                        parts = line.strip().split()
                        if len(parts) == 2:
                            time = float(parts[0])
                            count = int(parts[1])
                            accumulated_count += count
                            accumulated_time += time
                            times.append(accumulated_time)
                            counts.append(accumulated_count)

                    times.reverse()
                    counts.reverse()

                    if counts and times:
                        color = colors[i % len(colors)]
                        ax.plot(times, counts, marker='o', color=color, label=os.path.basename(selected_file))
                    else:
                        messagebox.showerror("Ошибка", f"Файл {selected_file} не содержит данных для анализа")

            ax.set_xlabel('Время выполнения (секунды)')
            ax.set_ylabel('Количество простых чисел')
            ax.set_title('Анализ времени выполнения')
            ax.legend()
            ax.grid(True)


            save_dir = os.path.dirname(self.selected_files[0])
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            save_path = os.path.join(save_dir, f"analysis_plot_{timestamp}.png")
            self.figure.savefig(save_path)
            messagebox.showinfo("Сохранение графика", f"График сохранен: {save_path}")

            plt.show()
        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка при чтении файлов: {e}")

    def reset(self):
        self.selected_files = []
        self.selected_files_label.config(text="")
        self.result_label.config(text="")
        if self.figure:
            plt.close(self.figure)
            self.figure = None

class StepPrimeApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Поиск простых чисел с шагом")
        self.root.geometry("500x700")
        self.root.resizable(width=False, height=False)

        self.prime_calculator = PrimeCalculator()

        self.label_start = tk.Label(self.root, text="Начальное число:")
        self.label_start.pack(pady=10)

        self.entry_start = tk.Entry(self.root)
        self.entry_start.pack(pady=10)

        self.label_end = tk.Label(self.root, text="Конечное число:")
        self.label_end.pack(pady=10)

        self.entry_end = tk.Entry(self.root)
        self.entry_end.pack(pady=10)

        self.label_step = tk.Label(self.root, text="Шаг:")
        self.label_step.pack(pady=10)

        self.entry_step = tk.Entry(self.root)
        self.entry_step.pack(pady=10)

        self.label_method = tk.Label(self.root, text="Выберите метод:")
        self.label_method.pack(pady=10)

        self.methods = [
            "Решето Эратосфена",
            "Решето Аткина",
            "Простой перебор",
            "Сегментированное решето",
            "Пробное деление",
            "Факторизация колеса"
        ]
        self.selected_method = tk.StringVar()
        self.selected_method.set(self.methods[0])

        for method in self.methods:
            tk.Radiobutton(self.root, text=method, variable=self.selected_method, value=method).pack()

        self.button_count = tk.Button(self.root, text="Подсчет", command=self.count_primes)
        self.button_count.pack(pady=10)

        self.button_analyze = tk.Button(self.root, text="Анализ результатов", command=self.analyze_results)
        self.button_analyze.pack(pady=10)

        self.button_exit = tk.Button(self.root, text="Выход", command=self.root.destroy)
        self.button_exit.pack(pady=10)

    def count_primes(self):
        try:
            start = int(self.entry_start.get())
            end = int(self.entry_end.get())
            step = int(self.entry_step.get())
        except ValueError:
            messagebox.showerror("Ошибка", "Пожалуйста, введите корректные числа")
            return

        if start >= end or step <= 0:
            messagebox.showerror("Ошибка", "Пожалуйста, убедитесь, что начальное число меньше конечного и шаг положителен")
            return

        method_selected = self.selected_method.get()

        primes = self.prime_calculator.find_primes_simple(start, end, step, method_selected)

        filename = f"result_primes_{method_selected}_{start}_{end}_{step}.txt"
        self.prime_calculator.write_to_file(primes, filename)

        result_str = f"Найдены простые числа в диапазоне от {start} до {end} с шагом {step} и методом {method_selected}\n"

        messagebox.showinfo("Результаты", result_str)

    def analyze_results(self):
        self.root.destroy()
        root_analysis = tk.Tk()
        app_analysis = AnalysisApp(root_analysis)
        root_analysis.mainloop()



root = tk.Tk()
app = PrimeApp(root)
root.mainloop()
