import { useEffect, useRef } from "react";
import type { ReactNode } from "react";

/** Native disclosure semantics, with outside-click and Escape dismissal. */
export function Popover({
  label,
  className = "",
  children,
  icon,
}: {
  label: string;
  className?: string;
  children: ReactNode;
  icon?: ReactNode;
}) {
  const ref = useRef<HTMLDetailsElement>(null);
  useEffect(() => {
    const dismiss = (event: PointerEvent) => {
      if (
        ref.current?.open &&
        event.target instanceof Node &&
        !ref.current.contains(event.target)
      )
        ref.current.open = false;
    };
    document.addEventListener("pointerdown", dismiss);
    return () => document.removeEventListener("pointerdown", dismiss);
  }, []);
  return (
    <details
      className={`popover ${className}`}
      ref={ref}
      onKeyDown={(event) => {
        if (event.key === "Escape" && ref.current?.open) {
          ref.current.open = false;
          ref.current.querySelector("summary")?.focus();
          event.stopPropagation();
        }
      }}
      onClick={(event) => {
        if (
          event.target instanceof Element &&
          event.target.closest("a[href]") &&
          ref.current
        )
          ref.current.open = false;
      }}
    >
      <summary aria-label={label} title={icon ? label : undefined}>
        {icon ?? label}
      </summary>
      <div className="popover-content">{children}</div>
    </details>
  );
}
