int main() {
#if __x86_64__
  return 1;
#else
  return 0;
#endif
}
