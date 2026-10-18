export function ErrorMessage({ error }: { error: Error | null }) {
  return error ? (
    <div className="error" role="alert">
      {error.message}
    </div>
  ) : null;
}
export function Loading({
  children = "Loading…",
}: {
  children?: React.ReactNode;
}) {
  return (
    <p className="muted" role="status">
      {children}
    </p>
  );
}
export const formatCount = (value: number) => value.toLocaleString();
